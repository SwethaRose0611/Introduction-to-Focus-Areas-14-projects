#include <divsufsort.h>
#include <sstream>
#include <stdexcept>
#include <chrono>

#include <seqan3/std/filesystem>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

// comparison between substrings
// saidx_t inherits from int32_t
int compare(std::vector<seqan3::dna5>& ref, size_t pos, size_t len, const std::vector<seqan3::dna5>& s)
{
    if (pos > ref.size())
    {
        throw std::out_of_range("Pos > Reference");
    } 

    for (size_t i = 0; i <= len; i++)
    {
        if (ref[i+pos] == s[i]) continue;
        if (ref[i+pos] <  s[i]) return -1;
        if (ref[i+pos] >  s[i]) return 1;
    }
    
    if (pos+len > ref.size()) return 1;
    if (pos+len < ref.size()) return -1;

    return 0;
}

void mlr_find(std::vector<seqan3::dna5>& query, std::vector<saidx_t>& sa, std::vector<seqan3::dna5>& ref)//, std::vector<uint32_t>& hits) 
{
    if (query.empty() || sa.empty() || ref.empty()) return;
    size_t pl = query.size(); //length of substring

    if (compare(ref, sa[0], pl, query) > 0) return; //query is smaller than smallest suffix
    if (compare(ref, sa[ref.size()-1], pl, query) < 0) return; //query is longer than longest suffix

    // left border
    size_t lp;  
    size_t m; size_t l = 0; size_t r = ref.size() - 1;

    while (r - l > 1) {
        m = (l + r + 1) / 2;    
        if (compare(ref, sa[m], pl, query) >= 0) 
            r = m; 
        else l = m;
    }
    if (compare(ref, sa[l], pl, query) != 0) 
        lp = r; 
    else lp = l;

    // right border
    size_t rp;
    l = 0; r = ref.size() - 1;

    while (r - l > 1) {
        m = (l + r + 1) / 2;    
        if (compare(ref, sa[m], pl, query) <= 0)
            l = m; 
        else r = m;
    }
    if (compare(ref, sa[r], pl, query) != 0) 
        rp = l; 
    else rp = r;

    // print hits
    for (size_t i = rp; i <= lp; i++) 
        std::cout << "\tFound query between " << sa[i] << " and " << sa[i]+query.size() << "\n";
}


int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"suffixarray_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    try {
         parser.parse();
    } catch (seqan3::argument_parser_error const& ext) {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // loading our files
    auto reference_stream = seqan3::sequence_file_input{reference_file};
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // read reference into memory
    // Attention: we are concatenating all sequences into one big combined sequence
    //            this is done to simplify the implementation of suffix_arrays
    std::vector<seqan3::dna5> reference;
    for (auto& record : reference_stream) {
        auto r = record.sequence();
        reference.insert(reference.end(), r.begin(), r.end());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    // Array that should hold the future suffix array
    std::vector<saidx_t> suffixarray;
    suffixarray.resize(reference.size()); // resizing the array, so it can hold the complete SA

    //!TODO !ImplementMe implement suffix array sort
    //Hint, if can use libdivsufsort (already integrated in this repo)
    //      https://github.com/y-256/libdivsufsort
    //      To make the `reference` compatible with libdivsufsort you can simply
    //      cast it by calling:
    //      `sauchar_t const* str = reinterpret_cast<sauchar_t const*>(reference.data());`
    sauchar_t const* str = reinterpret_cast<sauchar_t const*>(reference.data());
    divsufsort((unsigned char * )str, suffixarray.data(), reference.size());
    // output
    // for(int i = 0; i < suffixarray.size(); ++i) {
    //     printf("suffixarray[%2d] = %2d: ", i, suffixarray[i]);
    //     for(int j = suffixarray[i]; j < suffixarray.size(); ++j) {
    //         printf("%c", str[j]);
    //     }
    //     printf("$\n");
    // }

    // 100000 for task 5, since there was no difference to a lesser amount of queries
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    std::vector<size_t> query_numbers{100};//{1000, 10000, 100000, 1000000};
    for (size_t n_queries : query_numbers)
    {        
        //!TODO here adjust the number of searches
        queries.resize(n_queries); // will reduce the amount of searches
        // change to 1,000; 10,000; 100,000; 1,000,000

        // Start Search
        start = std::chrono::steady_clock::now();
        for (auto& q : queries) {
            //!TODO !ImplementMe apply binary search and find q  in reference using binary search on `suffixarray`
            // You can choose if you want to use binary search based on "naive approach", "mlr-trick", "lcp"
            // We used mlr
            mlr_find(q, suffixarray, reference);
        }
        // End of search and print time
        end = std::chrono::steady_clock::now();
        std::cout << "The Program took " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
            << "ms with " 
            << n_queries << " queries. \n";
    }
    return 0;
}
