// 원래 영역의 각 variant occurrence를 같은 read/HP/CIGAR 규칙으로 계수함.
#include <htslib/sam.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
struct Variant {
    int64_t begin,end,pos;std::string ref,alt,kind;
    uint64_t depth=0,alt_count=0,h1_alt=0,h1_depth=0,h2_alt=0,h2_depth=0;
};
struct Block {int64_t begin,end,query;};
struct Event {int op;int64_t pos,length;};
static int haplotype(bam1_t* b) {
    const uint8_t* tag=bam_aux_get(b,"HP");if(!tag)return 0;
    double value=0;
    if(std::string("cCsSiI").find(char(*tag))!=std::string::npos)value=bam_aux2i(tag);
    // 기존 float HP1/2의 dict key 오류를 조용히 정수로 바꾸지 않음.
    else if(*tag=='f'||*tag=='d'){value=bam_aux2f(tag);return value==1||value==2?3:0;}
    return value==1?1:value==2?2:0;
}
int main(int argc,char**argv) {
    samFile* in=nullptr;sam_hdr_t* header=nullptr;hts_idx_t* index=nullptr;hts_itr_t* iter=nullptr;bam1_t* read=nullptr;int rc=1;
    try {
        if(argc!=7)throw std::runtime_error("bam reference|- contig variants.tsv min_mapq decode_threads");
        int mapq=std::stoi(argv[5]),threads=std::stoi(argv[6]);
        if(threads<0)throw std::runtime_error("invalid decoder count");
        std::ifstream source(argv[4]);if(!source)throw std::runtime_error("variant input missing");
        std::vector<Variant> variants;std::string line;
        int64_t left=INT64_MAX,right=0;
        while(std::getline(source,line)) {
            std::istringstream stream(line);Variant v;std::string extra;
            if(!(stream>>v.begin>>v.end>>v.pos>>v.ref>>v.alt>>v.kind)||(stream>>extra)||v.begin<0||v.end<=v.begin||v.pos<0)throw std::runtime_error("invalid variant row");
            if(v.kind!="SNV"&&v.kind!="MNV"&&v.kind!="INDEL")throw std::runtime_error("invalid variant kind");
            if(v.kind!="INDEL"&&v.ref.size()!=v.alt.size())throw std::runtime_error("invalid MNV length");
            left=std::min(left,v.begin);right=std::max(right,v.end);variants.push_back(std::move(v));
        }
        if(!source.eof())throw std::runtime_error("variant input read error");
        std::vector<size_t> order(variants.size());std::iota(order.begin(),order.end(),0);
        std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return variants[a].pos<variants[b].pos;});
        in=sam_open(argv[1],"r");if(!in)throw std::runtime_error("alignment open failed");
        if(std::string(argv[2])!="-"&&hts_set_fai_filename(in,argv[2])!=0)throw std::runtime_error("reference setup failed");
        int eof=hts_check_EOF(in);if(eof==0||eof==-1)throw std::runtime_error("truncated alignment EOF");
        if(threads&&hts_set_threads(in,threads)!=0)throw std::runtime_error("decoder setup failed");
        header=sam_hdr_read(in);if(!header)throw std::runtime_error("alignment header failed");
        int tid=sam_hdr_name2tid(header,argv[3]);if(tid<0)throw std::runtime_error("contig missing");
        index=sam_index_load(in,argv[1]);if(!index)throw std::runtime_error("alignment index missing");
        uint64_t fetched=0,retained=0;
        if(!variants.empty()) {
            iter=sam_itr_queryi(index,tid,left,right);if(!iter)throw std::runtime_error("alignment iterator failed");
            read=bam_init1();if(!read)throw std::bad_alloc();int status;
            while((status=sam_itr_next(in,iter,read))>=0) {
                ++fetched;if((read->core.flag&3844)||int(read->core.qual)<mapq)continue;
                int64_t start=read->core.pos,finish=bam_endpos(read);
                auto first=std::lower_bound(order.begin(),order.end(),start-10,[&](size_t a,int64_t p){return variants[a].pos<p;});
                auto last=std::upper_bound(first,order.end(),finish,[&](int64_t p,size_t a){return p<variants[a].pos;});
                if(first==last)continue;++retained;
                std::vector<Block> blocks;std::vector<Event> events;int64_t ref=start,query=0;
                const uint32_t* cigar=bam_get_cigar(read);
                for(uint32_t i=0;i<read->core.n_cigar;++i) {
                    int op=bam_cigar_op(cigar[i]);int64_t n=bam_cigar_oplen(cigar[i]);
                    if(op==BAM_CBACK)throw std::runtime_error("unsupported backward CIGAR");
                    if(op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF)blocks.push_back({ref,ref+n,query});
                    if(op==BAM_CINS||op==BAM_CDEL||op==BAM_CREF_SKIP)events.push_back({op,ref,n});
                    if(bam_cigar_type(op)&2)ref+=n;
                    if(bam_cigar_type(op)&1)query+=n;
                    // 기존 pysam.get_aligned_pairs는 P도 query offset에 더함. 이 호환 동작을 유지함.
                    if(op==BAM_CPAD)query+=n;
                }
                int hp=haplotype(read);const uint8_t* sequence=bam_get_seq(read);
                for(auto p=first;p!=last;++p) {
                    Variant& v=variants[*p];
                    // VCF의 긴 DEL이 여러 원래 영역에 등장하는 경우도 별도로 보존함.
                    if(start>=v.end||finish<=v.begin)continue;
                    bool has=false,covered=true;
                    if(v.kind=="INDEL") {
                        if(start>v.pos-10||finish<v.pos+int64_t(v.ref.size())+10)continue;
                        int desired=v.ref.size()>v.alt.size()?BAM_CDEL:BAM_CINS;
                        int64_t size=std::llabs(int64_t(v.ref.size())-int64_t(v.alt.size()));
                        for(const Event& event:events) {
                            if(event.op==BAM_CREF_SKIP&&event.pos<=v.pos+1&&v.pos+1<event.pos+event.length){covered=false;break;}
                            if(event.op==desired&&event.length==size&&std::llabs(event.pos-(v.pos+1))<=10)has=true;
                        }
                    } else {
                        has=true;
                        for(size_t k=0;k<v.ref.size();++k) {
                            int64_t pos=v.pos+int64_t(k);
                            auto it=std::upper_bound(blocks.begin(),blocks.end(),pos,[](int64_t a,const Block& b){return a<b.begin;});
                            if(it==blocks.begin()){covered=false;break;}--it;
                            if(pos>=it->end){covered=false;break;}
                            int64_t q=it->query+pos-it->begin;
                            if(q<0||q>=read->core.l_qseq)throw std::runtime_error("query sequence/CIGAR mismatch");
                            char wanted=v.alt[k];if(wanted>='a'&&wanted<='z')wanted-=32;
                            if(seq_nt16_str[bam_seqi(sequence,q)]!=wanted)has=false;
                        }
                    }
                    if(!covered)continue;
                    if(hp==3)throw std::runtime_error("non-integer HP1/2 violates legacy counter keys");
                    ++v.depth;v.alt_count+=has;
                    if(hp==1){++v.h1_depth;v.h1_alt+=has;}else if(hp==2){++v.h2_depth;v.h2_alt+=has;}
                }
            }
            if(status<-1)throw std::runtime_error("truncated/corrupt alignment body");
        }
        for(const Variant& v:variants)std::cout<<v.depth<<'\t'<<v.alt_count<<'\t'<<v.h1_alt<<'\t'<<v.h1_depth<<'\t'<<v.h2_alt<<'\t'<<v.h2_depth<<'\n';
        if(!std::cout)throw std::runtime_error("count output failed");
        std::cerr<<"VARIANT_SUPPORT records="<<fetched<<" candidates="<<retained<<" variants="<<variants.size()<<'\n';rc=0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';}
    if(read)bam_destroy1(read);if(iter)hts_itr_destroy(iter);if(index)hts_idx_destroy(index);if(header)sam_hdr_destroy(header);if(in&&sam_close(in)<0)rc=1;return rc;
}
