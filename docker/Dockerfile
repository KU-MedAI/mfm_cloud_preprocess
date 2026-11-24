# syntax=docker/dockerfile:1
# Linux x86_64의 CPU 전처리 이미지임. 생물학 설정은 실행 JSON을 사용함.
FROM ubuntu:24.04@sha256:9c7f200c07e3d4a03e7249d2c16bf66a51120fb876aa2242cf31451412543289

SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG MINIFORGE_VERSION=26.7.2-0
ARG MINIFORGE_SHA256=281b0ac7d550802efc81af633225a5e6116d29ae72f3ab4eae7168c3931a4c05

# OS 업데이트는 허용하고 실제 배포 이미지는 digest로 고정함.
# hadolint ignore=DL3008
RUN test "$(uname -m)" = x86_64 \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       bash bzip2 ca-certificates curl libgomp1 procps tini util-linux \
    && rm -rf /var/lib/apt/lists/* \
    && curl --fail --location --retry 5 \
       "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh" \
       --output /tmp/miniforge.sh \
    && printf '%s  %s\n' "$MINIFORGE_SHA256" /tmp/miniforge.sh | sha256sum --check - \
    && bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh

# 환경 정의가 같으면 코드 수정 때에도 환경 설치 레이어를 재사용함.
COPY environment.yml requirements.txt /opt/mfm-env/
RUN /opt/conda/bin/conda env create --name jj_mfm --file /opt/mfm-env/environment.yml \
    && /opt/conda/envs/jj_mfm/bin/python -m pip install --no-cache-dir --no-deps -r /opt/mfm-env/requirements.txt \
    && /opt/conda/envs/jj_mfm/bin/python -m pip check \
    && /opt/conda/bin/conda clean --all --yes \
    && mkdir -p /opt/mfm-image \
    && dpkg-query -W > /opt/mfm-image/os-packages.tsv \
    && /opt/conda/bin/conda list --name jj_mfm --explicit > /opt/mfm-image/conda-explicit.txt \
    && /opt/conda/envs/jj_mfm/bin/python -m pip freeze > /opt/mfm-image/pip-freeze.txt

ENV MFM_CORE_ENV=/opt/conda/envs/jj_mfm \
    CONDA_PREFIX=/opt/conda/envs/jj_mfm \
    PATH=/opt/conda/envs/jj_mfm/bin:/opt/conda/bin:${PATH} \
    PYTHONPATH=/opt/mfm:/opt/mfm/src:/opt/mfm/vendor/MFM-DNARNAFM:/opt/mfm/tests \
    PYTHONDONTWRITEBYTECODE=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    HOME=/scratch/home \
    TMPDIR=/scratch/tmp \
    XDG_CACHE_HOME=/cache \
    MFM_NATIVE_CACHE_ROOT=/cache/native \
    MFM_IMAGE_NATIVE_CACHE=/opt/mfm-image/native

WORKDIR /opt/mfm
COPY . /opt/mfm/
RUN chmod -R a+rX /opt/mfm \
    && python tools/check_cloud_package.py \
    && bash 00_environment.sh verify-core \
    && bash 00_environment.sh verify-rna \
    && mkdir -p /opt/mfm-image/native /scratch/home /scratch/tmp /cache/native \
    && MFM_NATIVE_CACHE_ROOT=/opt/mfm-image/native python docker/smoke.py \
    && chmod -R a+rX /opt/mfm-image/native \
    && chmod 0755 docker/entrypoint.sh \
    && groupadd --gid 10001 mfm \
    && useradd --no-log-init --uid 10001 --gid 10001 --no-create-home --home-dir /scratch/home mfm \
    && chown -R mfm:mfm /scratch /cache \
    && mkdir -p /data/input /data/assets /data/output /data/state \
    && chown -R mfm:mfm /data

USER 10001:10001
STOPSIGNAL SIGTERM
ENTRYPOINT ["/usr/bin/tini", "-g", "--", "/opt/mfm/docker/entrypoint.sh"]
CMD ["help"]
