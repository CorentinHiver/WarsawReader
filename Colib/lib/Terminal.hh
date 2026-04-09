#pragma once

#include <iostream>

namespace Colib
{
  namespace Terminal
  {
    constexpr char CLEAR_ROW              [] = "\033[2K"  ;
    constexpr char DISABLE_LINES_WRAPPING [] = "\033[?7l" ;
    constexpr char ENABLE_LINES_WRAPPING  [] = "\033[?7h" ;
    constexpr char RETURN                 [] =        "\r";
    constexpr char NEWLINE                [] =        "\n";

    std::ostream& move_up   (int index = 1, std::ostream & out = std::cout) {out << "\033[" << index << "A"; return out;}
    std::ostream& move_down (int index = 1, std::ostream & out = std::cout) {out << "\033[" << index << "B"; return out;}
    
    std::ostream& goto_col (int colID, std::ostream & out = std::cout) {out << "\033[" << colID << "G"; return out;}

    std::ostream& clear_row (std::ostream & out = std::cout) {out << CLEAR_ROW; return out;}
    std::ostream& new_line  (std::ostream & out = std::cout) {out << NEWLINE  ; return out;}
    std::ostream& goto_left (std::ostream & out = std::cout) {out << RETURN   ; return out;}
    
    std::ostream& goto_right (std::ostream & out = std::cout) {out << "\033[" <<  999  << "G"; return out;}

    std::ostream& new_lines (int lines = 1, std::ostream & out = std::cout) {for (int i = 0; i<lines; ++i) out << NEWLINE; return out;}

    std::ostream& disable_line_wrapping (std::ostream & out = std::cout) {out << DISABLE_LINES_WRAPPING; return out;}
    std::ostream& enable_line_wrapping  (std::ostream & out = std::cout) {out << ENABLE_LINES_WRAPPING ; return out;}
    
    std::ostream& clear_row_return (std::ostream & out = std::cout) {out << CLEAR_ROW << RETURN; return out;}
    std::ostream& flush (std::ostream & out = std::cout) {out << std::flush; return out;}
  } 
}