// csBioLCCC.i - SWIG interface
%module csBioLCCC 
%include "std_string.i"

%{
#include "auxiliary.hpp"
#include "aminoacid.h"
#include "terminus.h"
#include "chemicalbasis.h"
#include "gradientpoint.h"
#include "gradient.h"
#include "chromoconditions.h"
#include "BioLCCC.h"
%}

// Parse the original header file
%include "auxiliary.hpp"
%include "aminoacid.h"
%include "terminus.h"
%include "chemicalbasis.h"
%include "gradientpoint.h"
%include "gradient.h"
%include "chromoconditions.h"
%include "BioLCCC.h"

// Instantiate some templates

