from __future__ import annotations

import hashlib
import json
import sqlite3
from pathlib import Path
from typing import Any


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def manifest_progress_signature(
    manifest_rows: list[dict[str, str]],
    config: dict[str, Any],
) -> str:
    """manifest 순서와 실행 설정이 모두 같은 경우에만 crash journal을 재사용한다."""
    payload = {
        "rows": manifest_rows,
        "config": _json_ready(config),
    }
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


class ManifestProgressJournal:
    """sample 완료를 O(1) SQLite transaction으로 기록한다.

    과거 구현은 sample 하나가 끝날 때마다 누적 status TSV와 전체 catalog를 다시
    썼기 때문에 N개 sample에서 O(N²) shared-filesystem I/O가 발생했다. journal은
    sample당 한 row만 즉시 commit하고 기존 TSV/catalog는 마지막에 한 번 materialize한다.
    비정상 종료 뒤 같은 manifest/config로 재실행하면 완료 sample을 복구한다.
    """

    def __init__(self, path: Path, signature: str) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(str(self.path), timeout=60.0)
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.execute("PRAGMA synchronous=FULL")
        self.db.execute(
            "CREATE TABLE IF NOT EXISTS meta "
            "(key TEXT PRIMARY KEY, value TEXT NOT NULL)"
        )
        self.db.execute(
            "CREATE TABLE IF NOT EXISTS result ("
            "seq INTEGER PRIMARY KEY AUTOINCREMENT,"
            "manifest_order INTEGER NOT NULL UNIQUE,"
            "sample_id TEXT NOT NULL,"
            "status TEXT NOT NULL,"
            "payload TEXT NOT NULL"
            ")"
        )
        row = self.db.execute(
            "SELECT value FROM meta WHERE key='signature'"
        ).fetchone()
        existing = str(row[0]) if row else ""
        if existing != signature:
            with self.db:
                self.db.execute("DELETE FROM result")
                self.db.execute(
                    "INSERT INTO meta(key,value) VALUES('signature',?) "
                    "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                    (signature,),
                )
        else:
            self.db.commit()

    def record(self, result: dict[str, Any]) -> None:
        payload = json.dumps(
            _json_ready(result),
            sort_keys=True,
            ensure_ascii=False,
            separators=(",", ":"),
        )
        with self.db:
            # 실패 row를 재시도한 경우 완료 시점의 completion order를 새로 부여한다.
            self.db.execute(
                "DELETE FROM result WHERE manifest_order=? AND status!='complete'",
                (int(result.get("manifest_order", 0)),),
            )
            self.db.execute(
                "INSERT INTO result(manifest_order,sample_id,status,payload) "
                "VALUES(?,?,?,?) "
                "ON CONFLICT(manifest_order) DO UPDATE SET "
                "sample_id=excluded.sample_id,status=excluded.status,"
                "payload=excluded.payload",
                (
                    int(result.get("manifest_order", 0)),
                    str(result.get("sample_id", "")),
                    str(result.get("status", "")),
                    payload,
                ),
            )

    def results(self, *, completion_order: bool = False) -> list[dict[str, Any]]:
        order = "seq" if completion_order else "manifest_order"
        rows = self.db.execute(
            f"SELECT payload FROM result ORDER BY {order}"
        ).fetchall()
        return [json.loads(str(row[0])) for row in rows]

    def completed_sample_ids(self) -> set[str]:
        return {
            str(row[0])
            for row in self.db.execute(
                "SELECT sample_id FROM result WHERE status='complete'"
            )
        }

    def close(self) -> None:
        if self.db is not None:
            self.db.close()
            self.db = None

    def cleanup(self) -> None:
        """정상 materialize 뒤 생성 journal만 제거한다."""
        if self.db is not None:
            self.db.execute("PRAGMA wal_checkpoint(TRUNCATE)")
            self.close()
        for path in (
            self.path,
            Path(str(self.path) + "-wal"),
            Path(str(self.path) + "-shm"),
        ):
            try:
                path.unlink()
            except FileNotFoundError:
                pass
