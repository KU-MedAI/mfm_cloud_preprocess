// 같은 QNAME의 mate/양HAP/교차염색체와 입력 순서를 같은 shard에 보존함.
#include <htslib/sam.h>
#include <cstdint>
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
template<class T>void put(std::ostream& out,T value){out.write(reinterpret_cast<const char*>(&value),sizeof(T));}
int main(int argc,char**argv) {
    samFile* in=nullptr;sam_hdr_t* h=nullptr;bam1_t* b=nullptr;int rc=1;
    try {
        if(argc<5 || argc>7)throw std::runtime_error("bam prefix shards duplicates");bool duplicates=std::stoi(argv[4])!=0;int minimum=argc==7?std::stoi(argv[6]):255;int shards=std::stoi(argv[3]);if(shards<1||shards>256)throw std::runtime_error("invalid shards");
        in=sam_open(argv[1],"r");if(!in)throw std::runtime_error("input missing");if(argc>=6 && std::stoi(argv[5])>0 && hts_set_threads(in,std::stoi(argv[5]))!=0)throw std::runtime_error("BGZF threads failed");int eof=hts_check_EOF(in);if(eof==0 || eof==-1)throw std::runtime_error("truncated BAM EOF");h=sam_hdr_read(in);if(!h)throw std::runtime_error("header failed");
        std::vector<std::unique_ptr<std::ofstream>> files;
        for(int i=0;i<shards;++i){auto out=std::make_unique<std::ofstream>(std::string(argv[2])+"."+std::to_string(i)+".bin",std::ios::binary);out->exceptions(std::ios::badbit|std::ios::failbit);out->write("MFMASE1\0",8);put(*out,uint32_t(sam_hdr_nref(h)));for(int t=0;t<sam_hdr_nref(h);++t){std::string name=sam_hdr_tid2name(h,t);put(*out,uint32_t(name.size()));out->write(name.data(),name.size());}files.push_back(std::move(out));}
        b=bam_init1();if(!b)throw std::bad_alloc();int status;uint64_t count=0;
        while((status=sam_read1(in,h,b))>=0) {
            if(b->core.flag&(BAM_FUNMAP|BAM_FSECONDARY|BAM_FSUPPLEMENTARY))continue;
            const char* name=bam_get_qname(b);std::string qname=name;uint64_t hash=14695981039346656037ull;for(unsigned char c:qname){hash^=c;hash*=1099511628211ull;}
            auto& out=*files[hash%shards];uint8_t* tag=((b->core.qual>0 || minimum<=0) && (duplicates || !(b->core.flag&BAM_FDUP)))?bam_aux_get(b,"AS"):nullptr;int64_t score=0;
            if(tag){if(*tag=='f'||*tag=='d'){double value=bam_aux2f(tag);if(!std::isfinite(value) || value>=double(INT64_MAX) || value<double(INT64_MIN))throw std::runtime_error("AS overflow");score=int64_t(value);}else if(std::string("cCsSiI").find(char(*tag))!=std::string::npos)score=bam_aux2i(tag);else if(*tag=='Z'||*tag=='A'){std::string value=*tag=='A'?std::string(1,bam_aux2A(tag)):std::string(bam_aux2Z(tag));size_t n;score=std::stoll(value,&n);if(value.find_first_not_of(" \t\r\n",n)!=std::string::npos)throw std::runtime_error("invalid AS text");}else throw std::runtime_error("invalid AS type");}
            std::vector<std::pair<int64_t,int64_t>> blocks;int64_t pos=b->core.pos;
            if(tag){const uint32_t* cigar=bam_get_cigar(b);for(uint32_t i=0;i<b->core.n_cigar;++i){int op=bam_cigar_op(cigar[i]);int64_t n=bam_cigar_oplen(cigar[i]);if(op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF)blocks.emplace_back(pos,pos+n);if(bam_cigar_type(op)&2)pos+=n;}}
            if(b->core.tid<0 || b->core.tid>=sam_hdr_nref(h))throw std::runtime_error("invalid mapped TID");
            put(out,uint32_t(qname.size()));put(out,uint16_t(b->core.flag));put(out,uint8_t(b->core.qual));put(out,uint8_t(tag!=nullptr));put(out,score);put(out,int32_t(b->core.tid));put(out,uint32_t(blocks.size()));out.write(qname.data(),qname.size());for(auto block:blocks){put(out,block.first);put(out,block.second);}++count;
        }
        if(status<-1)throw std::runtime_error("truncated/corrupt BAM");for(auto& f:files)f->close();std::cout<<"{\"records\":"<<count<<",\"shards\":"<<shards<<"}\n";rc=0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';}
    if(b)bam_destroy1(b);if(h)sam_hdr_destroy(h);if(in&&sam_close(in)<0)rc=1;return rc;
}
