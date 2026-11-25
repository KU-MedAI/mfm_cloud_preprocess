// BAM/CRAM을 한 번 읽어 flag·태그와 samtools 좌표 정렬 순서를 함께 확인함.
#include <htslib/sam.h>
#include <array>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <sys/resource.h>

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("alignment_scan input reference|- threads");
        int threads = std::stoi(argv[3]);
        if (threads < 0) throw std::runtime_error("음수 thread");
        auto started = std::chrono::steady_clock::now();
        std::unique_ptr<samFile, decltype(&hts_close)> in(sam_open(argv[1], "r"), hts_close);
        if (!in) throw std::runtime_error("입력 열기 실패");
        if (std::strcmp(argv[2], "-") && hts_set_fai_filename(in.get(), argv[2]) < 0)
            throw std::runtime_error("decode reference 설정 실패");
        if (threads && hts_set_threads(in.get(), threads) < 0)
            throw std::runtime_error("해독 thread 설정 실패");
        std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> header(sam_hdr_read(in.get()), sam_hdr_destroy);
        std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(bam_init1(), bam_destroy1);
        if (!header || !read) throw std::runtime_error("헤더/레코드 초기화 실패");
        kstring_t order = {0, 0, nullptr};
        bool declared = sam_hdr_find_tag_hd(header.get(), "SO", &order) == 0
                        && order.s && std::strcmp(order.s, "coordinate") == 0;
        std::free(order.s);
        std::array<uint64_t, 9> counts{};
        std::tuple<uint32_t, int64_t, bool> previous{};
        bool sorted = true, first = true;
        int status;
        while ((status = sam_read1(in.get(), header.get(), read.get())) >= 0) {
            const auto& core = read->core;
            // samtools 기본 sort: reference, position, reverse, 원래 입력 순서.
            auto key = std::make_tuple(uint32_t(core.tid), int64_t(core.pos), bool(core.flag & BAM_FREVERSE));
            if (!first && key < previous) sorted = false;
            previous = key; first = false;
            bool hp = false, ps = false;
            errno = 0;
            for (auto* aux = bam_aux_first(read.get()); aux; aux = bam_aux_next(read.get(), aux)) {
                const char* tag = bam_aux_tag(aux);
                hp |= tag[0] == 'H' && tag[1] == 'P';
                ps |= tag[0] == 'P' && tag[1] == 'S';
            }
            if (errno && errno != ENOENT) throw std::runtime_error("손상된 aux 태그");
            auto flag = core.flag;
            ++counts[0]; counts[1] += !(flag & BAM_FUNMAP);
            if (flag & BAM_FSECONDARY) ++counts[3];
            else if (flag & BAM_FSUPPLEMENTARY) ++counts[4];
            else ++counts[2];
            counts[5] += bool(flag & BAM_FDUP);
            counts[6] += bool(flag & BAM_FDUP) && !(flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));
            counts[7] += hp; counts[8] += ps;
        }
        if (status < -1) throw std::runtime_error("alignment 본문 해독 실패");
        if (hts_close(in.release()) != 0) throw std::runtime_error("alignment 닫기 실패");
        const char* names[] = {"records", "mapped_records", "primary_records", "secondary_records",
            "supplementary_records", "duplicate_records", "primary_duplicate_records",
            "hp_tagged_records", "ps_tagged_records"};
        struct rusage usage{}; getrusage(RUSAGE_SELF, &usage);
        double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        std::cout << "{\"counts\":{";
        for (size_t i = 0; i < counts.size(); ++i)
            std::cout << (i ? "," : "") << '"' << names[i] << "\":" << counts[i];
        std::cout << "},\"execution\":{\"passes\":1,\"native_threads\":" << 1 + int(threads > 0) + threads
                  << ",\"coordinate_sorted\":" << (sorted ? "true" : "false")
                  << ",\"header_coordinate\":" << (declared ? "true" : "false")
                  << ",\"seconds\":" << seconds << ",\"peak_rss_bytes\":" << uint64_t(usage.ru_maxrss) * 1024
                  << "}}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
