// 저장기의 CIGAR walker와 별개인 HTSlib pileup으로 reference 관측을 전수 집계함.
#include <htslib/sam.h>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

struct Output {
    uint32_t *bases;
    uint64_t *mapq, *baseq;
    uint32_t *baseq_n, *forward, *deletion, *insertion;
    uint32_t *deletion_event, *max_deletion, *supplementary;
};

struct Context {
    samFile *file = nullptr;
    sam_hdr_t *header = nullptr;
    hts_idx_t *index = nullptr;
    hts_itr_t *iterator = nullptr;
    int tid = -1, min_mapq = 0, excluded_flags = 0;
    bool proper_pair = false;
    int64_t start = 0, end = 0;
    Output *output = nullptr;
    std::string error;
    ~Context() {
        if (iterator) hts_itr_destroy(iterator);
        if (index) hts_idx_destroy(index);
        if (header) sam_hdr_destroy(header);
        if (file) sam_close(file);
    }
};

static void message(char *buffer, size_t size, const std::string &value) {
    if (size) std::snprintf(buffer, size, "%s", value.c_str());
}

static bool increment(uint32_t &value, Context *context) {
    if (value == UINT32_MAX) {
        context->error = "uint32 observation overflow";
        return false;
    }
    ++value;
    return true;
}

static int next_read(void *pointer, bam1_t *read) {
    auto *context = static_cast<Context *>(pointer);
    int status;
    while ((status = sam_itr_next(context->file, context->iterator, read)) >= 0) {
        const auto flags = read->core.flag;
        if ((flags & context->excluded_flags) || read->core.qual < context->min_mapq || !read->core.l_qseq)
            continue;
        if (context->proper_pair && (flags & BAM_FPAIRED) && !(flags & BAM_FPROPER_PAIR))
            continue;
        // I/D는 동일 anchor에 함께 올 수 있으므로 pileup의 단일 signed indel 값에 의존하지 않음.
        if (!(flags & BAM_FSUPPLEMENTARY)) {
            int64_t position = read->core.pos;
            const auto *cigar = bam_get_cigar(read);
            for (uint32_t i = 0; i < read->core.n_cigar;) {
                const int op = bam_cigar_op(cigar[i]);
                uint64_t length = 0;
                do { length += bam_cigar_oplen(cigar[i++]); }
                while (i < read->core.n_cigar && bam_cigar_op(cigar[i]) == op);
                if (op > BAM_CDIFF) {
                    context->error = "unsupported CIGAR operation";
                    return -2;
                }
                const int64_t anchor = position - 1;
                if ((op == BAM_CINS || op == BAM_CDEL) && position > read->core.pos &&
                    anchor >= context->start && anchor < context->end) {
                    auto &out = *context->output;
                    const auto offset = anchor - context->start;
                    if (!increment(op == BAM_CINS ? out.insertion[offset] : out.deletion_event[offset], context))
                        return -2;
                    if (op == BAM_CDEL) {
                        if (length > UINT32_MAX) { context->error = "deletion length overflow"; return -2; }
                        out.max_deletion[offset] = std::max(out.max_deletion[offset], static_cast<uint32_t>(length));
                    }
                }
                if (bam_cigar_type(op) & 2) position += length;
            }
        }
        return status;
    }
    if (status < -1) context->error = "alignment iterator read failure";
    return status;
}

extern "C" void *mfm_pileup_open(const char *path, const char *reference, const char *contig,
                                  int min_mapq, int excluded_flags, int proper_pair,
                                  char *error, size_t error_size) {
    try {
        auto context = std::make_unique<Context>();
        context->file = sam_open(path, "r");
        if (!context->file) throw std::runtime_error("alignment open failed");
        if (reference && *reference && hts_set_fai_filename(context->file, reference))
            throw std::runtime_error("reference setup failed");
        const int eof = hts_check_EOF(context->file);
        if (eof == 0 || eof == -1) throw std::runtime_error("alignment EOF check failed");
        context->header = sam_hdr_read(context->file);
        if (!context->header) throw std::runtime_error("alignment header failed");
        context->tid = sam_hdr_name2tid(context->header, contig);
        if (context->tid < 0) throw std::runtime_error("contig missing");
        context->index = sam_index_load(context->file, path);
        if (!context->index) throw std::runtime_error("alignment index missing");
        context->min_mapq = min_mapq;
        context->excluded_flags = excluded_flags;
        context->proper_pair = proper_pair;
        return context.release();
    } catch (const std::exception &exception) {
        message(error, error_size, exception.what());
        return nullptr;
    }
}

extern "C" int mfm_pileup_count(void *pointer, int64_t start, int64_t end, Output *out,
                                 char *error, size_t error_size) {
    auto *context = static_cast<Context *>(pointer);
    bam_plp_t pileup = nullptr;
    try {
        if (!context || !out || start < 0 || end <= start) throw std::runtime_error("invalid pileup range");
        context->start = start; context->end = end; context->output = out; context->error.clear();
        if (context->iterator) hts_itr_destroy(context->iterator);
        context->iterator = sam_itr_queryi(context->index, context->tid, start, end);
        if (!context->iterator) throw std::runtime_error("alignment iterator creation failed");
        pileup = bam_plp_init(next_read, context);
        if (!pileup) throw std::bad_alloc();
        bam_plp_set_maxcnt(pileup, INT_MAX);
        int tid, count;
        hts_pos_t position;
        const bam_pileup1_t *column;
        while ((column = bam_plp64_auto(pileup, &tid, &position, &count))) {
            if (position < start) continue;
            if (position >= end) break;
            if (count == INT_MAX) throw std::runtime_error("pileup integer depth limit reached; refusing silent cap");
            const auto offset = position - start;
            for (int i = 0; i < count; ++i) {
                const auto &entry = column[i];
                if (entry.is_refskip) continue;
                const auto *read = entry.b;
                if (read->core.flag & BAM_FSUPPLEMENTARY) {
                    if (!entry.is_del && !increment(out->supplementary[offset], context))
                        throw std::runtime_error(context->error);
                    continue;
                }
                if (entry.is_del) {
                    if (!increment(out->deletion[offset], context)) throw std::runtime_error(context->error);
                    continue;
                }
                if (entry.qpos < 0 || entry.qpos >= read->core.l_qseq) throw std::runtime_error("invalid query offset");
                const int base = bam_seqi(bam_get_seq(read), entry.qpos);
                const int code = base == 1 ? 0 : base == 2 ? 1 : base == 4 ? 2 : base == 8 ? 3 : base == 15 ? 4 : -1;
                if (code < 0) throw std::runtime_error("unsupported observed base");
                if (!increment(out->bases[code*(end-start)+offset], context)) throw std::runtime_error(context->error);
                out->mapq[offset] += read->core.qual;
                if (!(read->core.flag & BAM_FREVERSE) && !increment(out->forward[offset], context))
                    throw std::runtime_error(context->error);
                const auto *qualities = bam_get_qual(read);
                const uint8_t quality = qualities[entry.qpos];
                // pysam은 첫 QUAL이255이면 read.query_qualities 전체를 None으로 반환함.
                if (qualities[0] != 255 && quality != 255) {
                    out->baseq[offset] += quality;
                    if (!increment(out->baseq_n[offset], context)) throw std::runtime_error(context->error);
                }
            }
        }
        if (!context->error.empty()) throw std::runtime_error(context->error);
        bam_plp_destroy(pileup);
        return 0;
    } catch (const std::exception &exception) {
        if (pileup) bam_plp_destroy(pileup);
        message(error, error_size, exception.what());
        return 1;
    }
}

extern "C" void mfm_pileup_close(void *pointer) { delete static_cast<Context *>(pointer); }
