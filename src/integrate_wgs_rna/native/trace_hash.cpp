// canonical JSON 바이트를 보존하며 worker당 한 CPU로 배치 SHA256을 계산함.
#include <openssl/evp.h>
#include <stdexcept>
#include <charconv>
#include <cstdint>
#include <string>

struct Sha256 {
    EVP_MD_CTX* ctx=EVP_MD_CTX_new();
    Sha256() {if(!ctx)throw std::bad_alloc();}
    ~Sha256() {EVP_MD_CTX_free(ctx);}
    void digest(const std::string& text,uint8_t* out) {
        unsigned int size=0;
        if(EVP_DigestInit_ex(ctx,EVP_sha256(),nullptr)!=1 ||
           EVP_DigestUpdate(ctx,text.data(),text.size())!=1 ||
           EVP_DigestFinal_ex(ctx,out,&size)!=1 || size!=32)throw std::runtime_error("SHA256 failed");
    }
};

static void number(std::string& s, int64_t value) {
    char buf[32]; auto end=std::to_chars(buf,buf+sizeof(buf),value).ptr;
    s.append(buf,end);
}

extern "C" int trace_hash_batch(const uint8_t* bases, int64_t size,
    const int64_t* offsets, const int64_t* meta, const int64_t* cigar,
    int64_t cigar_size, int64_t count, int packed, uint8_t* out) {
    try {
        std::string json; Sha256 sha;
        for(int64_t i=0;i<count;++i) {
            const auto* m=meta+i*8;
            auto a=offsets[i], b=offsets[i+1], len=m[7], c=m[5], n=m[6];
            if(a<0 || b<a || b>size || len<0 || len>(b-a)*(packed?2:1)
                || c<0 || n<0 || c>cigar_size || n>cigar_size-c) return 1;
            json.clear(); json+="[\"";
            if(packed) {
                json.resize(2+len);
                for(int64_t j=0;j<len;++j) {
                    auto byte=bases[a+j/2];auto code=j%2?byte&15:byte>>4;
                    if(code>4)return 2;
                    json[2+j]="ACGTN"[code];
                }
            } else {
                int64_t from=0;
                for(int64_t j=0;j<len;++j) {
                    auto ch=bases[a+j];
                    if(ch<32 || ch>126)return 3;
                    if(ch=='"' || ch=='\\') {
                        json.append(reinterpret_cast<const char*>(bases+a+from),j-from);
                        json+='\\';json+=char(ch);from=j+1;
                    }
                }
                json.append(reinterpret_cast<const char*>(bases+a+from),len-from);
            }
            json+='"';
            for(int k=0;k<5;++k) {json+=',';number(json,m[k]);}
            json+=",[";
            for(int64_t j=0;j<n;++j) {if(j) json+=',';number(json,cigar[c+j]);}
            json+="]]";
            sha.digest(json,out+32*i);
        }
        return 0;
    } catch(...) {return 4;}
}

extern "C" int trace_segment_hash_batch(const uint8_t* reads, const int64_t* meta,
    int64_t count, uint8_t* out) {
    try {
        static const char hex[]="0123456789abcdef";
        std::string json; Sha256 sha;
        for(int64_t i=0;i<count;++i) {
            json="[\"";
            for(int j=0;j<32;++j) {auto ch=reads[i*32+j];json+=hex[ch>>4];json+=hex[ch&15];}
            json+='"';
            for(int j=0;j<6;++j) {json+=',';number(json,meta[i*6+j]);}
            json+=']';
            sha.digest(json,out+32*i);
        }
        return 0;
    } catch(...) {return 4;}
}
