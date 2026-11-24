// 인덱스 없이 BAM 전체를 한 번 읽고 완전한 레코드 묶음을 병렬 집계한다.
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <array>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <deque>
#include <exception>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <sys/resource.h>

using Bytes = std::vector<uint8_t>;
using Counts = std::array<uint64_t, 9>;
static uint32_t u32(const uint8_t* p) {
    return uint32_t(p[0]) | uint32_t(p[1]) << 8 | uint32_t(p[2]) << 16 | uint32_t(p[3]) << 24;
}
static void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }

static Counts count_batch(const Bytes& batch) {
    Counts counts{};
    size_t offset = 0;
    while (offset < batch.size()) {
        require(batch.size() - offset >= 4, "BAM 레코드 길이 누락");
        const size_t length = u32(batch.data() + offset);
        offset += 4;
        require(length >= 32 && length <= batch.size() - offset, "BAM 레코드 본문 손상");
        const uint8_t* p = batch.data() + offset;
        const size_t name = p[8], cigar = uint32_t(p[12]) | uint32_t(p[13]) << 8;
        const uint32_t sequence = u32(p + 16);
        require(sequence <= INT32_MAX && name > 0, "BAM sequence/name 길이 오류");
        const uint16_t flag = uint16_t(p[14]) | uint16_t(p[15]) << 8;
        size_t aux = 32 + name + 4 * cigar + (size_t(sequence) + 1) / 2 + sequence;
        require(aux <= length && p[32 + name - 1] == 0, "BAM sequence/quality 범위 오류");
        bool hp = false, ps = false;
        while (aux < length) {
            require(length - aux >= 3, "BAM aux 헤더 손상");
            const uint8_t a = p[aux], b = p[aux + 1], type = p[aux + 2];
            aux += 3;
            hp |= a == 'H' && b == 'P'; ps |= a == 'P' && b == 'S';
            size_t bytes = 0;
            switch (type) {
                case 'A': case 'c': case 'C': bytes = 1; break;
                case 's': case 'S': bytes = 2; break;
                case 'i': case 'I': case 'f': bytes = 4; break;
                case 'd': bytes = 8; break;
                case 'Z': case 'H': {
                    const void* end = std::memchr(p + aux, 0, length - aux);
                    require(end != nullptr, "BAM 문자열 aux 종료 누락");
                    bytes = static_cast<const uint8_t*>(end) - (p + aux) + 1;
                    break;
                }
                case 'B': {
                    require(length - aux >= 5, "BAM 배열 aux 헤더 손상");
                    size_t width = 0;
                    switch (p[aux]) {
                        case 'c': case 'C': width = 1; break;
                        case 's': case 'S': width = 2; break;
                        case 'i': case 'I': case 'f': width = 4; break;
                        default: throw std::runtime_error("BAM 배열 aux 타입 오류");
                    }
                    const uint32_t n = u32(p + aux + 1);
                    require(n <= INT32_MAX && n <= (length - aux - 5) / width, "BAM 배열 aux 길이 오류");
                    bytes = 5 + size_t(n) * width;
                    break;
                }
                default: throw std::runtime_error("BAM aux 타입 오류");
            }
            require(bytes <= length - aux, "BAM aux 본문 손상"); aux += bytes;
        }
        counts[0]++;
        counts[1] += !(flag & BAM_FUNMAP);
        if (flag & BAM_FSECONDARY) counts[3]++;
        else if (flag & BAM_FSUPPLEMENTARY) counts[4]++;
        else counts[2]++;
        counts[5] += bool(flag & BAM_FDUP);
        counts[6] += bool(flag & BAM_FDUP) && !(flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));
        counts[7] += hp; counts[8] += ps;
        offset += length;
    }
    return counts;
}

class Pipeline {
    std::mutex mutex;
    std::condition_variable changed;
    std::deque<Bytes> queue;
    std::vector<std::thread> workers;
    std::vector<Counts> totals;
    std::exception_ptr failure;
    bool done = false;
    size_t capacity;
public:
    explicit Pipeline(int count): totals(count ? count : 1), capacity(std::max(1, 2 * count)) {
        try {
            for (int i = 0; i < count; ++i) workers.emplace_back([this, i] {
                try {
                    while (true) {
                        Bytes batch;
                        {
                            std::unique_lock<std::mutex> lock(mutex);
                            changed.wait(lock, [this]{ return done || failure || !queue.empty(); });
                            if (failure || (done && queue.empty())) return;
                            batch = std::move(queue.front()); queue.pop_front(); changed.notify_all();
                        }
                        auto part = count_batch(batch);
                        for (size_t j = 0; j < part.size(); ++j) totals[i][j] += part[j];
                    }
                } catch (...) { abort(std::current_exception()); }
            });
        } catch (...) { abort(std::current_exception()); join(); throw; }
    }
    ~Pipeline() { abort(nullptr); join(); }
    void abort(std::exception_ptr error) {
        std::lock_guard<std::mutex> lock(mutex);
        if (error && !failure) failure = error;
        done = true; changed.notify_all();
    }
    void join() { for (auto& worker : workers) if (worker.joinable()) worker.join(); }
    void push(Bytes batch) {
        if (workers.empty()) {
            auto part = count_batch(batch);
            for (size_t i = 0; i < part.size(); ++i) totals[0][i] += part[i];
            return;
        }
        std::unique_lock<std::mutex> lock(mutex);
        changed.wait(lock, [this]{ return failure || queue.size() < capacity; });
        if (failure) std::rethrow_exception(failure);
        queue.push_back(std::move(batch)); changed.notify_all();
    }
    Counts finish() {
        { std::lock_guard<std::mutex> lock(mutex); done = true; changed.notify_all(); }
        join(); if (failure) std::rethrow_exception(failure);
        Counts result{};
        for (const auto& row : totals) for (size_t j = 0; j < row.size(); ++j) result[j] += row[j];
        return result;
    }
};

int main(int argc, char** argv) {
    BGZF* input = nullptr;
    try {
        require(argc == 6, "사용법: record_stats 입력 decode수 집계수 배치바이트 레코드상한");
        const int decode = std::stoi(argv[2]), workers = std::stoi(argv[3]);
        const size_t batch_size = std::stoull(argv[4]), record_limit = std::stoull(argv[5]);
        require(decode >= 0 && workers >= 0 && batch_size > 0 && record_limit >= 32, "자원 인자 오류");
        auto start = std::chrono::steady_clock::now();
        struct rusage initial{}; getrusage(RUSAGE_SELF, &initial);
        input = bgzf_open(argv[1], "r"); require(input != nullptr, "BAM 열기 실패");
        bam_hdr_t* header = bam_hdr_read(input); require(header != nullptr, "BAM 헤더 해독 실패"); bam_hdr_destroy(header);
        if (decode) require(bgzf_mt(input, decode, 64) == 0, "BGZF 병렬 해독 시작 실패");
        Pipeline pipeline(workers);
        size_t native_threads = 0;
        for (const auto& entry : std::filesystem::directory_iterator("/proc/self/task")) { (void)entry; native_threads++; }
        Bytes buffer; buffer.reserve(batch_size + 4);
        while (true) {
            const size_t prior = buffer.size(); buffer.resize(prior + batch_size);
            const ssize_t n = bgzf_read(input, buffer.data() + prior, batch_size);
            require(n >= 0, "BGZF 해독/CRC 검사 실패"); buffer.resize(prior + size_t(n));
            size_t complete = 0;
            while (buffer.size() - complete >= 4) {
                const size_t length = u32(buffer.data() + complete);
                require(length >= 32 && length <= INT32_MAX, "BAM 레코드 길이 손상");
                if (length > record_limit) throw std::length_error("레코드가 RAM 사전계산 상한을 초과함");
                if (buffer.size() - complete - 4 < length) break;
                complete += length + 4;
            }
            if (complete) {
                Bytes tail(buffer.begin() + complete, buffer.end()); buffer.resize(complete);
                pipeline.push(std::move(buffer)); buffer = std::move(tail);
            }
            if (!n) { require(buffer.empty(), "BAM 마지막 레코드가 잘림"); break; }
        }
        const auto counts = pipeline.finish();
        const int closed = bgzf_close(input); input = nullptr;
        require(closed == 0, "BGZF 닫기 실패");
        struct rusage usage{}; getrusage(RUSAGE_SELF, &usage);
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        const double cpu = usage.ru_utime.tv_sec - initial.ru_utime.tv_sec + (usage.ru_utime.tv_usec - initial.ru_utime.tv_usec) / 1e6
                         + usage.ru_stime.tv_sec - initial.ru_stime.tv_sec + (usage.ru_stime.tv_usec - initial.ru_stime.tv_usec) / 1e6;
        const char* names[] = {"records", "mapped_records", "primary_records", "secondary_records", "supplementary_records", "duplicate_records", "primary_duplicate_records", "hp_tagged_records", "ps_tagged_records"};
        std::cout << "{\"counts\":{";
        for (size_t i = 0; i < counts.size(); ++i) std::cout << (i ? "," : "") << '"' << names[i] << "\":" << counts[i];
        std::cout << "},\"execution\":{\"wall_seconds\":" << wall << ",\"cpu_seconds\":" << cpu << ",\"average_used_cores\":" << cpu / wall
                  << ",\"peak_rss_bytes\":" << usage.ru_maxrss * 1024L << ",\"native_threads\":" << native_threads
                  << ",\"decode_threads\":" << decode << ",\"count_workers\":" << workers << "}}\n";
        return 0;
    } catch (const std::length_error& error) {
        if (input) bgzf_close(input); std::cerr << error.what() << '\n'; return 75;
    } catch (const std::exception& error) {
        if (input) bgzf_close(input); std::cerr << error.what() << '\n'; return 1;
    }
}
