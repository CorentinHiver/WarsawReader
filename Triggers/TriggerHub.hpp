#pragma once

#ifdef TRIGGER

#ifndef STRINGIFY
#define STRINGIFY1(x) #x
#define STRINGIFY(x) STRINGIFY1(x)
#endif

#pragma message ("TRIGGER CHOSE TO BE " STRINGIFY(TRIGGER)) // This is a normal preprocessor message

// See TriggerExample.hpp for examples

#include HEADER(TRIGGER) // This line triggers an error in VSCode but compiles just fine

#endif //TRIGGER