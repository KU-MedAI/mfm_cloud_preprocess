// gzip/plain FASTQ를 한 번 읽고 QNAME/SEQ/QUAL 보존·통계를 함께 계산함.
#include <zlib.h>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
static bool line(gzFile in,std::string& value) {
    value.clear();char block[65536];
    while(gzgets(in,block,sizeof(block))) {
        value+=block;if(value.size()>(64u<<20))throw std::runtime_error("FASTQ line too large");
        if(value.back()=='\n')return true;
    }
    int err=0;gzerror(in,&err);if(err!=Z_OK && err!=Z_STREAM_END)throw std::runtime_error("corrupt gzip input");
    return !value.empty();
}
static size_t bases(const std::string& value) {size_t n=value.size();while(n && (value[n-1]=='\n'||value[n-1]=='\r'))--n;return n;}
int main(int argc,char**argv) {
    gzFile in=nullptr;FILE* out=nullptr;int rc=1;
    try {
        if(argc!=4)throw std::runtime_error("input destination|- mate(0/1/2)");
        in=gzopen(argv[1],"rb");if(!in)throw std::runtime_error("FASTQ input missing");gzbuffer(in,1<<20);
        if(std::strcmp(argv[2],"-")){out=fopen(argv[2],"wb");if(!out)throw std::runtime_error("output open failed");setvbuf(out,nullptr,_IOFBF,1<<20);}
        int mate=std::stoi(argv[3]);if(mate<0 || mate>2)throw std::runtime_error("invalid mate");
        uint64_t count=0,maxlen=0;std::string name,seq,plus,qual;
        while(line(in,name)) {
            if(!line(in,seq)||!line(in,plus)||!line(in,qual)||name[0]!='@'||plus[0]!='+'||bases(seq)!=bases(qual))throw std::runtime_error("invalid/truncated FASTQ");
            name.resize(bases(name));
            if(mate){std::string suffix="/"+std::to_string(mate);if(name.size()<2 || name.compare(name.size()-2,2,suffix))throw std::runtime_error("mate suffix missing");name.resize(name.size()-2);}
            name+='\n';
            if(out)for(const auto* value:{&name,&seq,&plus,&qual})if(fwrite(value->data(),1,value->size(),out)!=value->size())throw std::runtime_error("FASTQ write failed");
            // 기존 awk는 LF만 제거하므로 CRLF의 CR도 길이에 포함함.
            uint64_t len=seq.size()-(seq.back()=='\n');if(len>maxlen)maxlen=len;++count;
        }
        if(gzclose(in)!=Z_OK){in=nullptr;throw std::runtime_error("gzip close failed");}in=nullptr;
        if(out){if(fclose(out)){out=nullptr;throw std::runtime_error("output close failed");}out=nullptr;}
        std::cout<<"{\"records\":"<<count<<",\"max_read_length\":"<<maxlen<<",\"backend\":\"native\"}\n";rc=0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';}
    if(in)gzclose(in);if(out)fclose(out);return rc;
}
