#include "ParameterInput.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

std::string ParameterInput::trim(const std::string &str) {
  const auto first = std::find_if_not(str.begin(), str.end(),
      [](unsigned char c) {
        return std::isspace(c);
      });

  if (first == str.end()) {
    return "";
  }

  const auto last = std::find_if_not(str.rbegin(), str.rend(),
      [](unsigned char c) {
        return std::isspace(c);
      }).base();

  return std::string(first, last);
}


void ParameterInput::addFromString(const std::string& input) {
    std::istringstream stream(input);
    std::string token;
    std::string currentHeading;

    while (std::getline(stream, token, ',')) {
        token = trim(token);

        if (token.empty()) {
            continue;
        }

        const auto pos = token.find('=');

        // A token without '=' is a heading.
        if (pos == std::string::npos) {
            currentHeading = token;
            continue;
        }

        // Ignore key-value pairs before the first heading.
        if (currentHeading.empty()) {
            continue;
        }

        const std::string key = trim(token.substr(0, pos));
        const std::string value = trim(token.substr(pos + 1));

        if (!key.empty()) {
            blocks_[currentHeading][key] = value;
        }
    }
}



bool ParameterInput::load(const std::string &filename) {
  std::ifstream file(filename);

  if (!file.is_open()) {
    return false;
  }

  PIBlocks loaded;
  std::string currentHeading;
  std::string line;

  while (std::getline(file, line)) {
    line = trim(line);

    // Ignore blank lines and comments.
    if (line.empty() || line[0] == '#' || line[0] == ';') {
      continue;
    }

    // Detect headings: <Heading>
    if (line.front() == '<' && line.back() == '>') {

      currentHeading = trim(line.substr(1, line.size() - 2));

      if (!currentHeading.empty()) {
        loaded.try_emplace(currentHeading);
      } else {
        currentHeading.clear();
      }

      continue;
    }

    // Ignore parameters outside a heading.
    if (currentHeading.empty()) {
      continue;
    }

    // Parse key=value.
    const auto pos = line.find('=');

    if (pos == std::string::npos) {
      continue;
    }

    std::string key = trim(line.substr(0, pos));
    std::string value = trim(line.substr(pos + 1));

    if (!key.empty()) {
      loaded[currentHeading][key] = value;
    }
  }

  if (file.bad()) {
    return false;
  }

  blocks_ = std::move(loaded);
  return true;
}

std::string ParameterInput::get(const std::string &heading,
    const std::string &key) const {

  const auto blockIt = blocks_.find(heading);

  if (blockIt == blocks_.end()) {
    throw std::out_of_range(
        "ParameterInput: heading '" + heading + "' does not exist.");
  }

  const auto paramIt = blockIt->second.find(key);

  if (paramIt == blockIt->second.end()) {
    throw std::out_of_range(
        "ParameterInput: key '" + key + "' does not exist in heading '"
            + heading + "'.");
  }

  return paramIt->second;
}

std::string ParameterInput::getAndAdd(const std::string &heading,
    const std::string &key, const std::string &addValue) {

  auto &parameters = blocks_[heading];

  const auto it = parameters.find(key);

  if (it != parameters.end()) {
    return it->second;
  }

  return parameters.emplace(key, addValue).first->second;
}

bool ParameterInput::contains(const std::string &heading,
    const std::string &key) const {

  const auto blockIt = blocks_.find(heading);

  return blockIt != blocks_.end()
      && blockIt->second.find(key) != blockIt->second.end();
}
Real ParameterInput::getReal(const std::string &heading,
    const std::string &key) const {

  return static_cast<Real>(std::stold(get(heading, key).c_str()));
}

Real ParameterInput::getAndAddReal(const std::string &heading,
    const std::string &key, Real addValue) {

  const std::string value = getAndAdd(heading, key, std::to_string(addValue));

  return static_cast<Real>(std::stold(value.c_str()));
}
const ParameterInput::PIBlocks& ParameterInput::values() const {
  return blocks_;
}
