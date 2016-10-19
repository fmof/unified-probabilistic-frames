#ifndef FERRUM_UNUSED_VAR_HPP_
#define FERRUM_UNUSED_VAR_HPP_

// from http://stackoverflow.com/questions/12198449/cross-platform-macro-for-silencing-unused-variables-warning/12199209#12199209
#define MON_Internal_UnusedStringify(macro_arg_string_literal) #macro_arg_string_literal

#define MONUnusedParameter(macro_arg_parameter) _Pragma(MON_Internal_UnusedStringify(unused(macro_arg_parameter)))

template<typename T>
inline void ignore_result(T /* unused result */) {
}

#endif
