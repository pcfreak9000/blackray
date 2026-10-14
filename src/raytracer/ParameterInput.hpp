#pragma once

#include <map>
#include <string>

#include "def.hpp"

class ParameterInput {
public:
    using PIParameters = std::map<std::string, std::string>;
    using PIBlocks = std::map<std::string, PIParameters>;

    bool load(const std::string& filename);

    void addFromString(const std::string& input);

    std::string get(
        const std::string& heading,
        const std::string& key) const;

    std::string getAndAdd(
        const std::string& heading,
        const std::string& key,
        const std::string& addValue);

    // Real-valued overloads
    Real getReal(
        const std::string& heading,
        const std::string& key) const;

    Real getAndAddReal(
        const std::string& heading,
        const std::string& key,
        Real addValue);

    bool contains(
        const std::string& heading,
        const std::string& key) const;

    const PIBlocks& values() const;

private:
    PIBlocks blocks_;

    static std::string trim(const std::string& str);
};
