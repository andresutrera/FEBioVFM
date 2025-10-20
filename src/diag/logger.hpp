// #pragma once
// #include <string>
// enum class LogLevel { Error, Warn, Info, Debug, Trace };

// struct LogSink {
//   // Default: FEBio sink
//   static void set_level(LogLevel lvl);
//   static LogLevel level();
//   static void write(LogLevel lvl, const char* fmt, ...); // routes to feLog*/feLogDebugEx
// };

// // convenience
// #define VFM_ERROR(fmt, ...) ::LogSink::write(LogLevel::Error, fmt, ##__VA_ARGS__)
// #define VFM_WARN(fmt, ...)  ::LogSink::write(LogLevel::Warn,  fmt, ##__VA_ARGS__)
// #define VFM_INFO(fmt, ...)  ::LogSink::write(LogLevel::Info,  fmt, ##__VA_ARGS__)
// #define VFM_DEBUG(fmt, ...) ::LogSink::write(LogLevel::Debug, fmt, ##__VA_ARGS__)
// #define VFM_TRACE(fmt, ...) ::LogSink::write(LogLevel::Trace, fmt, ##__VA_ARGS__)