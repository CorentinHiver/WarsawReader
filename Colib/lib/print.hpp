#pragma once

// This is used in my PC to have better looking errors
#ifndef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0/1
#endif //_GLIBCXX_USE_CXX11_ABI

#include <iostream>
#include <iomanip>
#include <vector>

#include "Terminal.hh"

/// Stream vector support
template <class T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) noexcept
{
  for (const auto & e : v) os << e << " ";
  return os;
}

namespace Colib
{
  namespace Color
  {
    constexpr char RESET         [] = "\u001b[0m";
  
    constexpr char BLACK         [] = "\u001b[30m";
    constexpr char RED           [] = "\u001b[31m";
    constexpr char GREEN         [] = "\u001b[32m";
    constexpr char YELLOW        [] = "\u001b[33m";
    constexpr char BLUE          [] = "\u001b[34m";
    constexpr char MAGENTA       [] = "\u001b[35m";
    constexpr char CYAN          [] = "\u001b[36m";
    constexpr char WHITE         [] = "\u001b[37m";
  
    constexpr char BRIGHTBLACK   [] = "\u001b[30;1m";
    constexpr char BRIGHTRED     [] = "\u001b[31;1m";
    constexpr char BRIGHTGREEN   [] = "\u001b[32;1m";
    constexpr char BRIGHTYELLOW  [] = "\u001b[33;1m";
    constexpr char BRIGHTBLUE    [] = "\u001b[34;1m";
    constexpr char BRIGHTMAGENTA [] = "\u001b[35;1m";
    constexpr char BRIGHTCYAN    [] = "\u001b[36;1m";
    constexpr char BRIGHTWHITE   [] = "\u001b[37;1m";
  
    constexpr char GREY          [] = "\u001b[38;5;8m";
  }
}

#ifdef CoMT
inline std::mutex print_mutex;
#endif // CoMT

// using MT::lock_mutex = std::lock_guard<std::mutex>;

/// Print newline only
inline void print() noexcept
{
#ifdef CoMT
  Colib::MT::inject_print([](){ std::cout << Colib::Terminal::NEWLINE; });
#else // no CoMT
  std::cout << std::endl;
#endif // CoMT
}

/// Print with space separation + newline
template <class... Args>
constexpr inline void print(const Args&... args) noexcept
{
#ifdef CoMT
  Colib::MT::inject_print([&]() {
    ((std::cout << args << ' '), ...);
    std::cout << Colib::Terminal::NEWLINE;
  });
#else // no CoMT
  ((std::cout << args << ' '), ...);
  std::cout << std::endl;
#endif // CoMT
}

/// @brief Concatenated print + newline
template <class... Args>
constexpr inline void printc(const Args&... args) noexcept
{
#ifdef CoMT
  Colib::MT::inject_print([&]() {((std::cout << args), ...); Colib::Terminal::new_line();});
#else // no CoMT
    ((std::cout << args), ...); std::cout << std::endl;
#endif // CoMT
}

/// @brief Print with space separation, no newline
template <class... Args>
constexpr inline void println(const Args&... args) noexcept
{
#ifdef CoMT
  Colib::MT::inject_print([&]() {((std::cout << args << ' '), ...);});
#else // no CoMT
  ((std::cout << args << ' '), ...);
  std::cout << std::flush;
#endif // CoMT

}

/// @brief Concatenated print, no newline - equivalent to std::cout
template <class... Args>
constexpr inline void printcln(const Args&... args) noexcept
{
#ifdef CoMT
  Colib::MT::inject_print([&]() {((std::cout << args), ...);});
#else // no CoMT
  ((std::cout << args), ...);
#endif // CoMT
}

/// @brief Print in the same line. In MT mode, returns true everytime
template <class... Ts> 
constexpr inline void printsln(Ts &&... ts)  noexcept
{
#ifdef CoMT
  Colib::MT::printsln(std::forward<Ts>(ts)...);
#else // no CoMT
  Colib::Terminal::clear_row();
  Colib::Terminal::goto_left();
  println(std::forward<Ts>(ts)...); 
#endif // CoMT
}

/// Debug print (print only if -DDEBUG is passed in the compilation line)
template <class... Args>
constexpr inline void debug([[maybe_unused]] Args &&... args) noexcept
{
#ifdef DEBUG
  print(std::forward<Args>(args)...);
#endif
}

/// @brief Print in red
template <class... T> constexpr inline void error      (T const & ... t) noexcept { print(Colib::Color::RED   , t..., Colib::Color::RESET); }
/// @brief Print in yellow
template <class... T> constexpr inline void warning    (T const & ... t) noexcept { print(Colib::Color::YELLOW, t..., Colib::Color::RESET); }
/// @brief Print in grey
template <class... T> constexpr inline void information(T const & ... t) noexcept {print(Colib::Color::GREY   , t..., Colib::Color::RESET); }

#ifdef Cpp20
constexpr inline std::string nicer_bool(bool const some_bool) noexcept
#else // Cpp<20
inline std::string nicer_bool(bool const some_bool) noexcept
#endif // Cpp20
{
  using namespace Colib::Color;
  if (some_bool) return (BLUE + std::string("true ") + RESET);
  else           return (RED  + std::string("false") + RESET);
}

/// @brief Prints a fixed-size output filled with 0 where t is not printed
template <size_t __length__, char fillWith = '0', class... T>
inline constexpr void printfill(T const &... t) noexcept {
  (printcln(std::setfill(fillWith), std::setw(__length__), t, " "), ...);
  print();
}

// Hextech :

template <class... T>
inline void printh(T &&... t) {print(std::hex, t..., std::dec);}

template <size_t length, char fillWith = '0', class... T>
inline void printhfill(T &&... t)
{
  printcln(std::hex);
  printfill<length, fillWith>(t...);
  printc(std::dec);
}
