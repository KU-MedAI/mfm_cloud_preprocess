// RNA MAPQ255/primary/QC-pass, DUP 포함, D 포함/N 제외, mate 겹침 두 번임.
#include <htslib/sam.h>
#include <cstdint>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
int main(int argc,char**argv) {
    samFile* in=nullptr;sam_hdr_t* h=nullptr;hts_idx_t* idx=nullptr;hts_itr_t* it=nullptr;bam1_t* b=nullptr;int rc=1;
    try {
        if(argc!=5)throw std::runtime_error("bam contig length depth|difference");
        int64_t length=std::stoll(argv[3]);if(length<1)throw std::runtime_error("invalid length");
        in=sam_open(argv[1],"r");if(!in)throw std::runtime_error("input open failed");
        h=sam_hdr_read(in);if(!h)throw std::runtime_error("header failed");
        int tid=sam_hdr_name2tid(h,argv[2]);if(tid<0 || sam_hdr_tid2len(h,tid)!=length)throw std::runtime_error("contig/length mismatch");
        idx=sam_index_load(in,argv[1]);if(!idx)throw std::runtime_error("input index missing");
        it=sam_itr_queryi(idx,tid,0,length);if(!it)throw std::runtime_error("iterator failed");
        std::vector<int64_t> delta(length+1,0);b=bam_init1();if(!b)throw std::bad_alloc();
        int status;uint64_t fetched=0,retained=0;
        while((status=sam_itr_next(in,it,b))>=0) {
            ++fetched;if((b->core.flag&2820) || b->core.qual<255)continue;
            ++retained;int64_t pos=b->core.pos;const uint32_t* cigar=bam_get_cigar(b);
            for(uint32_t j=0;j<b->core.n_cigar;++j) {
                int op=bam_cigar_op(cigar[j]);int64_t size=bam_cigar_oplen(cigar[j]);
                if(op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF || op==BAM_CDEL) {
                    int64_t lo=std::max<int64_t>(0,pos),hi=std::min(length,pos+size);
                    if(hi>lo){++delta[lo];--delta[hi];}
                }
                if(bam_cigar_type(op)&2)pos+=size;
            }
        }
        if(status<-1)throw std::runtime_error("truncated/corrupt alignment");
        std::string mode=argv[4];
        if(mode=="difference") {
            if(fwrite(delta.data(),sizeof(int64_t),delta.size(),stdout)!=delta.size())throw std::runtime_error("difference write failed");
        } else if(mode=="depth") {
            std::vector<uint32_t> output(1<<20);size_t used=0;int64_t value=0;
            for(int64_t i=0;i<length;++i) {
                value+=delta[i];if(value<0 || uint64_t(value)>UINT32_MAX)throw std::runtime_error("depth uint32 overflow");
                output[used++]=uint32_t(value);
                if(used==output.size()){if(fwrite(output.data(),4,used,stdout)!=used)throw std::runtime_error("depth write failed");used=0;}
            }
            if(used && fwrite(output.data(),4,used,stdout)!=used)throw std::runtime_error("depth write failed");
            if(value+delta[length]!=0)throw std::runtime_error("unbalanced depth");
        } else throw std::runtime_error("invalid output mode");
        if(fflush(stdout))throw std::runtime_error("output close failed");
        std::cerr<<"RNA_DEPTH fetched="<<fetched<<" retained="<<retained<<'\n';rc=0;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';}
    if(b)bam_destroy1(b);if(it)hts_itr_destroy(it);if(idx)hts_idx_destroy(idx);if(h)sam_hdr_destroy(h);if(in && sam_close(in)<0)rc=1;
    return rc;
}
