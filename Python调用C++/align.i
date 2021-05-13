%module align

%include "cpointer.i"
%include "carrays.i"

%pointer_class(double,double_p);
%array_class(double,double_array);

%include"align.h"

%{
#include"align.h"

%}

