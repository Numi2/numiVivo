#include "NumiVivoCore/Core.hpp"

#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>

namespace nvivo::json {

Value::Value() noexcept : storage_(nullptr) {}
Value::Value(std::nullptr_t) noexcept : storage_(nullptr) {}
Value::Value(bool value) noexcept : storage_(value) {}
Value::Value(double value) noexcept : storage_(value) {}
Value::Value(std::string value) : storage_(std::move(value)) {}
Value::Value(Array value) : storage_(std::move(value)) {}
Value::Value(Object value) : storage_(std::move(value)) {}

bool Value::isNull() const noexcept { return std::holds_alternative<std::nullptr_t>(storage_); }
bool Value::isBool() const noexcept { return std::holds_alternative<bool>(storage_); }
bool Value::isNumber() const noexcept { return std::holds_alternative<double>(storage_); }
bool Value::isString() const noexcept { return std::holds_alternative<std::string>(storage_); }
bool Value::isArray() const noexcept { return std::holds_alternative<Array>(storage_); }
bool Value::isObject() const noexcept { return std::holds_alternative<Object>(storage_); }

bool Value::asBool(bool fallback) const noexcept {
    if (const auto* value = std::get_if<bool>(&storage_)) {
        return *value;
    }
    return fallback;
}

double Value::asNumber(double fallback) const noexcept {
    if (const auto* value = std::get_if<double>(&storage_)) {
        return *value;
    }
    return fallback;
}

std::string_view Value::asString() const noexcept {
    if (const auto* value = std::get_if<std::string>(&storage_)) {
        return *value;
    }
    return {};
}

const Value::Array& Value::asArray() const {
    static const Array empty;
    if (const auto* value = std::get_if<Array>(&storage_)) {
        return *value;
    }
    return empty;
}

const Value::Object& Value::asObject() const {
    static const Object empty;
    if (const auto* value = std::get_if<Object>(&storage_)) {
        return *value;
    }
    return empty;
}

const Value* Value::get(std::string_view key) const noexcept {
    const auto* object = std::get_if<Object>(&storage_);
    if (object == nullptr) {
        return nullptr;
    }
    const auto iterator = object->find(key);
    return iterator == object->end() ? nullptr : &iterator->second;
}

const Value* Value::at(std::size_t index) const noexcept {
    const auto* array = std::get_if<Array>(&storage_);
    if (array == nullptr || index >= array->size()) {
        return nullptr;
    }
    return &(*array)[index];
}

const Value::Storage& Value::storage() const noexcept { return storage_; }

namespace {

class Parser {
public:
    Parser(std::string_view input, ParseLimits limits)
        : input_(input), limits_(limits) {}

    ParseResult run() {
        ParseResult result;
        if (input_.size() > limits_.maxBytes) {
            diagnostics_.fatal(
                "NVJ001",
                "JSON document exceeds the configured byte limit.",
                "$",
                "Reduce the source document or raise the explicit parser limit."
            );
            result.diagnostics = std::move(diagnostics_);
            return result;
        }

        skipWhitespace();
        auto root = parseValue(0);
        skipWhitespace();

        if (root.has_value() && position_ != input_.size()) {
            addError("NVJ002", "Unexpected trailing content after the JSON document.");
            root.reset();
        }

        if (!diagnostics_.hasErrors()) {
            result.root = std::move(root);
        }
        result.diagnostics = std::move(diagnostics_);
        return result;
    }

private:
    std::optional<Value> parseValue(std::size_t depth) {
        if (depth > limits_.maxDepth) {
            addError("NVJ003", "JSON nesting exceeds the configured depth limit.");
            return std::nullopt;
        }
        if (++nodeCount_ > limits_.maxNodes) {
            addError("NVJ004", "JSON node count exceeds the configured limit.");
            return std::nullopt;
        }
        if (position_ >= input_.size()) {
            addError("NVJ005", "Unexpected end of JSON input.");
            return std::nullopt;
        }

        switch (input_[position_]) {
            case '{': return parseObject(depth + 1);
            case '[': return parseArray(depth + 1);
            case '"': {
                auto value = parseString();
                if (!value.has_value()) {
                    return std::nullopt;
                }
                return Value(std::move(*value));
            }
            case 't': return parseLiteral("true", Value(true));
            case 'f': return parseLiteral("false", Value(false));
            case 'n': return parseLiteral("null", Value(nullptr));
            default:
                if (input_[position_] == '-' || isDigit(input_[position_])) {
                    return parseNumber();
                }
                addError("NVJ006", "Unexpected token while parsing a JSON value.");
                return std::nullopt;
        }
    }

    std::optional<Value> parseObject(std::size_t depth) {
        advance(); // {
        skipWhitespace();
        Value::Object object;

        if (consume('}')) {
            return Value(std::move(object));
        }

        while (position_ < input_.size()) {
            if (object.size() >= limits_.maxObjectMembers) {
                addError("NVJ007", "JSON object exceeds the configured member limit.");
                return std::nullopt;
            }
            if (input_[position_] != '"') {
                addError("NVJ008", "Expected a quoted object member name.");
                return std::nullopt;
            }

            auto key = parseString();
            if (!key.has_value()) {
                return std::nullopt;
            }
            skipWhitespace();
            if (!consume(':')) {
                addError("NVJ009", "Expected ':' after an object member name.");
                return std::nullopt;
            }
            skipWhitespace();

            auto value = parseValue(depth);
            if (!value.has_value()) {
                return std::nullopt;
            }

            const auto [iterator, inserted] = object.emplace(std::move(*key), std::move(*value));
            if (!inserted) {
                addError("NVJ010", "Duplicate object member names are not permitted.");
                return std::nullopt;
            }

            skipWhitespace();
            if (consume('}')) {
                return Value(std::move(object));
            }
            if (!consume(',')) {
                addError("NVJ011", "Expected ',' or '}' in a JSON object.");
                return std::nullopt;
            }
            skipWhitespace();
        }

        addError("NVJ012", "Unterminated JSON object.");
        return std::nullopt;
    }

    std::optional<Value> parseArray(std::size_t depth) {
        advance(); // [
        skipWhitespace();
        Value::Array array;

        if (consume(']')) {
            return Value(std::move(array));
        }

        while (position_ < input_.size()) {
            if (array.size() >= limits_.maxArrayElements) {
                addError("NVJ013", "JSON array exceeds the configured element limit.");
                return std::nullopt;
            }

            auto value = parseValue(depth);
            if (!value.has_value()) {
                return std::nullopt;
            }
            array.push_back(std::move(*value));

            skipWhitespace();
            if (consume(']')) {
                return Value(std::move(array));
            }
            if (!consume(',')) {
                addError("NVJ014", "Expected ',' or ']' in a JSON array.");
                return std::nullopt;
            }
            skipWhitespace();
        }

        addError("NVJ015", "Unterminated JSON array.");
        return std::nullopt;
    }

    std::optional<std::string> parseString() {
        if (!consume('"')) {
            addError("NVJ016", "Expected a JSON string.");
            return std::nullopt;
        }

        std::string output;
        while (position_ < input_.size()) {
            const char character = input_[position_];
            advance();

            if (character == '"') {
                if (output.size() > limits_.maxStringBytes) {
                    addError("NVJ017", "JSON string exceeds the configured byte limit.");
                    return std::nullopt;
                }
                return output;
            }

            if (static_cast<unsigned char>(character) < 0x20U) {
                addError("NVJ018", "Unescaped control character in JSON string.");
                return std::nullopt;
            }

            if (character != '\\') {
                output.push_back(character);
                if (output.size() > limits_.maxStringBytes) {
                    addError("NVJ017", "JSON string exceeds the configured byte limit.");
                    return std::nullopt;
                }
                continue;
            }

            if (position_ >= input_.size()) {
                addError("NVJ019", "Unterminated escape sequence in JSON string.");
                return std::nullopt;
            }

            const char escapeCode = input_[position_];
            advance();
            switch (escapeCode) {
                case '"': output.push_back('"'); break;
                case '\\': output.push_back('\\'); break;
                case '/': output.push_back('/'); break;
                case 'b': output.push_back('\b'); break;
                case 'f': output.push_back('\f'); break;
                case 'n': output.push_back('\n'); break;
                case 'r': output.push_back('\r'); break;
                case 't': output.push_back('\t'); break;
                case 'u': {
                    auto scalar = parseHexQuad();
                    if (!scalar.has_value()) {
                        return std::nullopt;
                    }

                    std::uint32_t codePoint = *scalar;
                    if (codePoint >= 0xD800U && codePoint <= 0xDBFFU) {
                        if (position_ + 2 > input_.size() || input_[position_] != '\\' ||
                            input_[position_ + 1] != 'u') {
                            addError("NVJ020", "High UTF-16 surrogate is not followed by a low surrogate.");
                            return std::nullopt;
                        }
                        advance();
                        advance();
                        auto low = parseHexQuad();
                        if (!low.has_value() || *low < 0xDC00U || *low > 0xDFFFU) {
                            addError("NVJ021", "Invalid UTF-16 low surrogate in JSON string.");
                            return std::nullopt;
                        }
                        codePoint = 0x10000U + ((codePoint - 0xD800U) << 10U) + (*low - 0xDC00U);
                    } else if (codePoint >= 0xDC00U && codePoint <= 0xDFFFU) {
                        addError("NVJ022", "Unexpected UTF-16 low surrogate in JSON string.");
                        return std::nullopt;
                    }

                    appendUTF8(output, codePoint);
                    break;
                }
                default:
                    addError("NVJ023", "Invalid JSON string escape sequence.");
                    return std::nullopt;
            }
        }

        addError("NVJ024", "Unterminated JSON string.");
        return std::nullopt;
    }

    std::optional<std::uint32_t> parseHexQuad() {
        if (position_ + 4 > input_.size()) {
            addError("NVJ025", "Incomplete Unicode escape sequence.");
            return std::nullopt;
        }

        std::uint32_t value = 0;
        for (int index = 0; index < 4; ++index) {
            const char character = input_[position_];
            advance();
            value <<= 4U;
            if (character >= '0' && character <= '9') {
                value |= static_cast<std::uint32_t>(character - '0');
            } else if (character >= 'a' && character <= 'f') {
                value |= static_cast<std::uint32_t>(character - 'a' + 10);
            } else if (character >= 'A' && character <= 'F') {
                value |= static_cast<std::uint32_t>(character - 'A' + 10);
            } else {
                addError("NVJ026", "Invalid hexadecimal digit in Unicode escape sequence.");
                return std::nullopt;
            }
        }
        return value;
    }

    std::optional<Value> parseNumber() {
        const std::size_t start = position_;

        consume('-');
        if (position_ >= input_.size()) {
            addError("NVJ027", "Incomplete JSON number.");
            return std::nullopt;
        }

        if (consume('0')) {
            if (position_ < input_.size() && isDigit(input_[position_])) {
                addError("NVJ028", "Leading zero is not permitted in a JSON number.");
                return std::nullopt;
            }
        } else {
            if (!isDigitOneToNine(input_[position_])) {
                addError("NVJ029", "Invalid integer component in JSON number.");
                return std::nullopt;
            }
            while (position_ < input_.size() && isDigit(input_[position_])) {
                advance();
            }
        }

        if (consume('.')) {
            if (position_ >= input_.size() || !isDigit(input_[position_])) {
                addError("NVJ030", "Fractional JSON number requires at least one digit.");
                return std::nullopt;
            }
            while (position_ < input_.size() && isDigit(input_[position_])) {
                advance();
            }
        }

        if (position_ < input_.size() && (input_[position_] == 'e' || input_[position_] == 'E')) {
            advance();
            if (position_ < input_.size() && (input_[position_] == '+' || input_[position_] == '-')) {
                advance();
            }
            if (position_ >= input_.size() || !isDigit(input_[position_])) {
                addError("NVJ031", "Exponent in JSON number requires at least one digit.");
                return std::nullopt;
            }
            while (position_ < input_.size() && isDigit(input_[position_])) {
                advance();
            }
        }

        const std::string token(input_.substr(start, position_ - start));
        errno = 0;
        char* end = nullptr;
        const double value = std::strtod(token.c_str(), &end);
        if (errno == ERANGE || end != token.c_str() + token.size() || !std::isfinite(value)) {
            addError("NVJ032", "JSON number is outside the supported finite range.");
            return std::nullopt;
        }
        return Value(value);
    }

    std::optional<Value> parseLiteral(std::string_view expected, Value value) {
        if (input_.substr(position_, expected.size()) != expected) {
            addError("NVJ033", "Invalid JSON literal.");
            return std::nullopt;
        }
        for (std::size_t index = 0; index < expected.size(); ++index) {
            advance();
        }
        return value;
    }

    void skipWhitespace() {
        while (position_ < input_.size()) {
            const char character = input_[position_];
            if (character != ' ' && character != '\t' && character != '\n' && character != '\r') {
                break;
            }
            advance();
        }
    }

    bool consume(char expected) {
        if (position_ < input_.size() && input_[position_] == expected) {
            advance();
            return true;
        }
        return false;
    }

    void advance() {
        if (position_ >= input_.size()) {
            return;
        }
        if (input_[position_] == '\n') {
            ++line_;
            column_ = 1;
        } else {
            ++column_;
        }
        ++position_;
    }

    void addError(std::string code, std::string message) {
        diagnostics_.add({
            Severity::error,
            std::move(code),
            std::move(message),
            "$",
            {},
            {position_, 1, line_, column_}
        });
    }

    static bool isDigit(char character) noexcept {
        return character >= '0' && character <= '9';
    }

    static bool isDigitOneToNine(char character) noexcept {
        return character >= '1' && character <= '9';
    }

    static void appendUTF8(std::string& output, std::uint32_t codePoint) {
        if (codePoint <= 0x7FU) {
            output.push_back(static_cast<char>(codePoint));
        } else if (codePoint <= 0x7FFU) {
            output.push_back(static_cast<char>(0xC0U | (codePoint >> 6U)));
            output.push_back(static_cast<char>(0x80U | (codePoint & 0x3FU)));
        } else if (codePoint <= 0xFFFFU) {
            output.push_back(static_cast<char>(0xE0U | (codePoint >> 12U)));
            output.push_back(static_cast<char>(0x80U | ((codePoint >> 6U) & 0x3FU)));
            output.push_back(static_cast<char>(0x80U | (codePoint & 0x3FU)));
        } else {
            output.push_back(static_cast<char>(0xF0U | (codePoint >> 18U)));
            output.push_back(static_cast<char>(0x80U | ((codePoint >> 12U) & 0x3FU)));
            output.push_back(static_cast<char>(0x80U | ((codePoint >> 6U) & 0x3FU)));
            output.push_back(static_cast<char>(0x80U | (codePoint & 0x3FU)));
        }
    }

    std::string_view input_;
    ParseLimits limits_;
    Diagnostics diagnostics_;
    std::size_t position_ = 0;
    std::size_t line_ = 1;
    std::size_t column_ = 1;
    std::size_t nodeCount_ = 0;
};

} // namespace

ParseResult parse(std::string_view input, const ParseLimits& limits) {
    return Parser(input, limits).run();
}

std::string escape(std::string_view value) {
    static constexpr char hexadecimal[] = "0123456789abcdef";
    std::string output;
    output.reserve(value.size() + 8);

    for (const unsigned char character : value) {
        switch (character) {
            case '"': output += "\\\""; break;
            case '\\': output += "\\\\"; break;
            case '\b': output += "\\b"; break;
            case '\f': output += "\\f"; break;
            case '\n': output += "\\n"; break;
            case '\r': output += "\\r"; break;
            case '\t': output += "\\t"; break;
            default:
                if (character < 0x20U) {
                    output += "\\u00";
                    output.push_back(hexadecimal[(character >> 4U) & 0xFU]);
                    output.push_back(hexadecimal[character & 0xFU]);
                } else {
                    output.push_back(static_cast<char>(character));
                }
        }
    }
    return output;
}

} // namespace nvivo::json
