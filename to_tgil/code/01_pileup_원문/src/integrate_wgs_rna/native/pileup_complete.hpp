// 공통 수식은 Python에서 계산하고 native는 관측만 집계함.
#pragma once
#include <htslib/sam.h>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

namespace complete_pileup {
enum Field { A, C, G, T, N, BQ_N, FORWARD, DELETION, INSERTION, DELETION_EVENT,
    MAX_DELETION, SUPPLEMENTARY, FIELDS };
static_assert(FIELDS==12,"count ABI must match schema.COUNT_FIELDS");
struct Output { uint32_t *counts; uint64_t *sums; };
inline void increment(uint32_t &value) {
    if (value==UINT32_MAX) throw std::overflow_error("uint32 observation overflow");
    ++value;
}
inline void add(uint64_t &value, uint64_t amount) {
    if (UINT64_MAX-value<amount) throw std::overflow_error("uint64 quality sum overflow");
    value+=amount;
}
struct Context {
    samFile *file=nullptr; sam_hdr_t *header=nullptr; hts_idx_t *index=nullptr; hts_itr_t *iterator=nullptr;
    int tid=-1, min_mapq=0, excluded_flags=0; bool proper_pair=false;
    int64_t start=0,end=0; Output output{}; std::string error;
    ~Context() {
        if(iterator)hts_itr_destroy(iterator);
        if(index)hts_idx_destroy(index);
        if(header)sam_hdr_destroy(header);
        if(file)sam_close(file);
    }
    uint32_t &count(int group,int field,int64_t position) {
        return output.counts[(int64_t(group)*FIELDS+field)*(end-start)+position-start];
    }
    uint64_t &sum(int group,int field,int64_t position) {
        return output.sums[(int64_t(group)*2+field)*(end-start)+position-start];
    }
    bool accepted(const bam1_t *read) const {
        int flag=read->core.flag;
        return !(flag&excluded_flags) && read->core.qual>=min_mapq && read->core.l_qseq>0 &&
               !(proper_pair && (flag&BAM_FPAIRED) && !(flag&BAM_FPROPER_PAIR));
    }
    static int hp(const bam1_t *read) {
        uint8_t *tag=bam_aux_get(read,"HP");
        if(!tag || std::string("cCsSiI").find(char(*tag))==std::string::npos)return 0;
        int64_t value=bam_aux2i(tag);return value==1 || value==2 ? int(value):0;
    }
    template<class Function> void groups(int phase,Function function) {
        function(0);if(phase)function(phase);
    }
    void base(const bam1_t *read,int64_t position,int64_t query,int phase) {
        if(position<start || position>=end)return;
        if(query<0 || query>=read->core.l_qseq)throw std::runtime_error("invalid query offset");
        bool supplementary=read->core.flag&BAM_FSUPPLEMENTARY;
        int code=bam_seqi(bam_get_seq(read),query),field=-1;
        if(!supplementary) {
            switch(code){case 1:field=A;break;case 2:field=C;break;case 4:field=G;break;case 8:field=T;break;case 15:field=N;break;}
            if(field<0)throw std::runtime_error("unsupported observed read base");
        }
        uint8_t quality=bam_get_qual(read)[query];
        groups(phase,[&](int group){
            if(supplementary){increment(count(group,SUPPLEMENTARY,position));return;}
            increment(count(group,field,position));
            if(!(read->core.flag&BAM_FREVERSE))increment(count(group,FORWARD,position));
            if(group==0) {
                add(sum(group,0,position),read->core.qual);
                if(quality!=255){increment(count(group,BQ_N,position));add(sum(group,1,position),quality);}
            }
        });
    }
    void deletion(const bam1_t *read,int64_t position,int phase) {
        if(position<start || position>=end || read->core.flag&BAM_FSUPPLEMENTARY)return;
        groups(phase,[&](int group){increment(count(group,DELETION,position));});
    }
    void event(const bam1_t *read,int64_t anchor,int64_t length,bool insertion,int64_t query,int phase) {
        if(anchor<start || anchor>=end || read->core.flag&BAM_FSUPPLEMENTARY)return;
        if(length<=0 || uint64_t(length)>UINT32_MAX)throw std::runtime_error("invalid indel length");
        if(insertion && (query<0 || query+length>read->core.l_qseq))throw std::runtime_error("insertion exceeds query");
        groups(phase,[&](int group){
            increment(count(group,insertion ? INSERTION:DELETION_EVENT,anchor));
            if(!insertion)count(group,MAX_DELETION,anchor)=std::max(count(group,MAX_DELETION,anchor),uint32_t(length));
        });
    }
    void cigar(const bam1_t *read,bool events_only) {
        int phase=hp(read);
        int64_t ref=read->core.pos,query=0;const uint32_t *cigar=bam_get_cigar(read);
        for(uint32_t i=0;i<read->core.n_cigar;) {
            int op=bam_cigar_op(cigar[i]);uint64_t length=0;
            do{length+=bam_cigar_oplen(cigar[i++]);}while(i<read->core.n_cigar && bam_cigar_op(cigar[i])==op);
            if(!length || op>BAM_CDIFF)throw std::runtime_error("invalid CIGAR operation");
            if(op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF) {
                if(query+int64_t(length)>read->core.l_qseq)throw std::runtime_error("CIGAR exceeds query");
                if(!events_only)for(int64_t p=std::max(ref,start);p<std::min(ref+int64_t(length),end);++p)base(read,p,query+p-ref,phase);
                query+=length;ref+=length;
            } else if(op==BAM_CINS) {
                if(ref>read->core.pos)event(read,ref-1,length,true,query,phase);
                query+=length;
            } else if(op==BAM_CDEL) {
                if(ref>read->core.pos)event(read,ref-1,length,false,query,phase);
                if(!events_only)for(int64_t p=std::max(ref,start);p<std::min(ref+int64_t(length),end);++p)deletion(read,p,phase);
                ref+=length;
            } else if(op==BAM_CREF_SKIP)ref+=length;
            else if(op==BAM_CSOFT_CLIP)query+=length;
        }
        if(query!=read->core.l_qseq)throw std::runtime_error("CIGAR query length mismatch");
    }
    void range(int64_t a,int64_t b,Output *out) {
        if(!out || a<0 || b<=a)throw std::runtime_error("invalid count range");
        start=a;end=b;output=*out;error.clear();
        if(iterator){hts_itr_destroy(iterator);iterator=nullptr;}
        iterator=sam_itr_queryi(index,tid,start,end);
        if(!iterator)throw std::runtime_error("alignment iterator creation failed");
    }
};
inline void message(char *out,size_t size,const std::string &value){if(size)std::snprintf(out,size,"%s",value.c_str());}
inline Context *open(const char *path,const char *reference,const char *contig,int mq,int flags,int proper) {
    auto context=std::make_unique<Context>();context->file=sam_open(path,"r");
    if(!context->file)throw std::runtime_error("alignment open failed");
    if(reference && *reference && hts_set_fai_filename(context->file,reference))throw std::runtime_error("reference setup failed");
    int eof=hts_check_EOF(context->file);if(eof==0 || eof==-1)throw std::runtime_error("alignment EOF check failed");
    context->header=sam_hdr_read(context->file);if(!context->header)throw std::runtime_error("alignment header failed");
    context->tid=sam_hdr_name2tid(context->header,contig);if(context->tid<0)throw std::runtime_error("contig missing");
    context->index=sam_index_load(context->file,path);if(!context->index)throw std::runtime_error("alignment index missing");
    context->min_mapq=mq;context->excluded_flags=flags;context->proper_pair=proper;return context.release();
}
} // namespace complete_pileup
