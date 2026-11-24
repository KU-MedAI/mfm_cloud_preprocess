// FASTQ 바이트/태그를 보존하고 검증된 mate별 최대 4MiB 프레임을 출력한다.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

constexpr size_t LIMIT = 64u << 20, BATCH = 4u << 20;
struct Input {
    std::array<char, 1u << 20> data;
    size_t pos = 0, end = 0;
    bool line(std::string &out) {
        out.clear();
        while (true) {
            if (pos == end) {
                end = fread(data.data(), 1, data.size(), stdin); pos = 0;
                if (!end) {
                    if (ferror(stdin)) throw std::runtime_error("FASTQ 입력 읽기 실패");
                    return !out.empty();
                }
            }
            auto newline = static_cast<const char *>(memchr(data.data()+pos, '\n', end-pos));
            size_t size = newline ? size_t(newline-data.data()-pos)+1 : end-pos;
            if (out.size()+size > LIMIT) throw std::length_error("FASTQ 단일 행이 bounded 입력 버퍼를 초과함");
            out.append(data.data()+pos, size); pos += size;
            if (newline) return true;
        }
    }
};
void write_bytes(const void *data, size_t size) {
    if (fwrite(data, 1, size, stdout) != size) throw std::runtime_error("FASTQ 프레임 출력 실패");
}
void frame(unsigned char mate, const char *data, uint32_t size) {
    unsigned char header[5] = {mate};
    for (int i=0; i<4; ++i) header[i+1] = (size >> (8*i)) & 255;
    write_bytes(header, 5); write_bytes(data, size);
}
size_t trimmed(const std::string &s) {
    size_t n=s.size(); while(n && (s[n-1]=='\n' || s[n-1]=='\r')) --n; return n;
}
std::string name(const std::string &s) {
    std::string n=s.substr(0, s.find_first_of(" \t\r\n\v\f"));
    for (const char *suffix : {"/1", "/2"})
        if (n.size()>=2 && n.compare(n.size()-2,2,suffix)==0) n.resize(n.size()-2);
    return n;
}
int main(int argc, char **argv) {
    try {
        if (argc!=2 || (std::string(argv[1])!="0" && std::string(argv[1])!="1")) throw std::runtime_error("paired 인자가 필요함");
        bool paired=std::string(argv[1])=="1";
        Input in; std::array<std::string,4> record;
        std::vector<std::string> buffers(paired ? 2 : 1);
        for (auto &buffer:buffers) buffer.reserve(BATCH);
        uint64_t count=0; std::string previous;
        while(in.line(record[0])) {
            for(int i=1;i<4;++i) in.line(record[i]);
            if(record[0][0]!='@' || record[2].empty() || record[2][0]!='+' || record[3].empty() || trimmed(record[1])!=trimmed(record[3]))
                throw std::invalid_argument("GPU 입력 FASTQ 레코드 손상");
            if(paired) {
                auto current=name(record[0]);
                if(count%2==0) previous=current;
                else if(current!=previous) throw std::invalid_argument("GPU paired FASTQ R1/R2 QNAME 불일치");
            }
            unsigned char mate=paired ? count%2 : 0;
            auto &buffer=buffers[mate];
            for(const auto &line:record) {
                size_t pos=0;
                while(pos<line.size()) {
                    size_t n=std::min(BATCH-buffer.size(), line.size()-pos);
                    buffer.append(line.data()+pos,n); pos+=n;
                    if(buffer.size()==BATCH) { frame(mate,buffer.data(),buffer.size()); buffer.clear(); }
                }
            }
            ++count;
        }
        if(paired && count%2) throw std::invalid_argument("GPU paired FASTQ mate 누락");
        if(!count) throw std::invalid_argument("GPU 입력 FASTQ가 비어 있음");
        for(size_t i=0;i<buffers.size();++i) if(!buffers[i].empty()) frame(i,buffers[i].data(),buffers[i].size());
        unsigned char encoded[8]; for(int i=0;i<8;++i) encoded[i]=(count>>(8*i))&255;
        frame(2,reinterpret_cast<const char *>(encoded),8);
        if(fflush(stdout)) throw std::runtime_error("FASTQ 프레임 flush 실패");
        return 0;
    } catch(const std::length_error &e) { fprintf(stderr,"%s\n",e.what()); return 75;
    } catch(const std::invalid_argument &e) { fprintf(stderr,"%s\n",e.what()); return 65;
    } catch(const std::exception &e) { fprintf(stderr,"%s\n",e.what()); return 1; }
}
