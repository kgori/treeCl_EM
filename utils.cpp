//
// Created by Kevin Gori on 15/06/15.
//

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "utils.h"

namespace utils{
    bool is_file(std::string filename) {
        std::ifstream fl(filename.c_str());
        bool result = true;
        if (!fl) {
            result = false;
        }
        fl.close();
        return result;
    }

    std::vector<std::string> readlines(const std::string& filename) {
        std::ifstream infile(filename);
        std::vector<std::string> result;
        std::string line;
        if (!infile) {
            std::cerr << "Error opening file" << std::endl;
            return result;
        }

        while (std::getline(infile, line)) result.push_back(line);
        return result;
    }

    alignmentUPtr parse_alignment_file(std::string path) {
        if (!is_file(path)) {
            std::cerr << "Couldn't find the alignment file " << path << std::endl;
            throw std::exception();
        }
        alignmentUPtr alignment;
        AlignmentDeleter del;
        alignment = alignmentUPtr(pllParseAlignmentFile(PLL_FORMAT_PHYLIP, path.c_str()), del);
        if (!alignment) {
            //std::cout << "Trying to parse as fasta" << std::endl;
            alignment = alignmentUPtr(pllParseAlignmentFile(PLL_FORMAT_FASTA, path.c_str()), del);
        }
        if (!alignment) {
            std::cerr << "Couldn't parse the alignment at " << path << std::endl;
            throw std::exception();
        }
        return alignment;
    }

    queueUPtr parse_partitions(std::string partitions) {
        if (is_file(partitions)) {
            return queueUPtr(pllPartitionParse(partitions.c_str()), QueueDeleter());
        }
        else {
            return queueUPtr(pllPartitionParseString(partitions.c_str()), QueueDeleter());
        }
    }

    double logsumexp(const std::vector<double>& nums) {
        double max_exp = nums[0], sum = 0.0;
        size_t i;

        for (i = 1 ; i < nums.size() ; i++)
            if (nums[i] > max_exp)
                max_exp = nums[i];

        for (i = 0; i < nums.size() ; i++)
            sum += exp(nums[i] - max_exp);

        return log(sum) + max_exp;
    }

    template<typename T>
    std::vector<T> csum(const std::vector<T> &row) {
        std::vector<T> result;
        double val = 0;
        for (double n : row) {
            val += n;
            result.push_back(val);
        }
        return result;
    }

    template std::vector<double> csum(const std::vector<double> &row);
    template std::vector<int> csum(const std::vector<int> &row);

    std::vector<double> scale_by_sum(const std::vector<double> &nums) {
        double sum = std::accumulate(nums.begin(), nums.end(), 0.0);
        std::cout << sum << std::endl;
        std::vector<double> output;
        for (auto elem : nums) {
            output.push_back(elem / sum);
        }
        return output;
    }

    size_t random_select(const std::vector<double> &probs) {
        std::random_device rd;
        std::mt19937 engine(rd());
        std::uniform_real_distribution<double> dist(0, 1);
        double u = dist(engine);
        for (size_t s = 0; s < probs.size(); ++s) {
            if (u < probs[s]) return s;
        }
        return probs.size();
    }
}
