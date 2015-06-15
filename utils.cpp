//
// Created by Kevin Gori on 15/06/15.
//

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
}
