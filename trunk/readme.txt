This library implements

- An N-dimensional SURF descriptor. SURF is somewhat similar to Sift, but can be computed faster because it makes use of rectangle filters in combination with integral images

- N-dimensional integral images.

- N-dimensional Gaussians and N-dimensional PCA.

See the paper [1] for more details on N-SURF.

It requires ITK.

Anyone is welcome to make additions and changes. Please contact 
jo <at> feulner-caps dot de 
if you need SVN write access.


License: MIT
#################
Copyright (C) 2012 Johannes Feulner

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
#################

Apart from this, it is cite-ware: Resulting scientific publications should cite [1].


Johannes Feulner, January 2012


[1] Comput Med Imaging Graph. 2011 Apr;35(3):227-36. 2010 Dec 3.
"Comparing axial CT slices in quantized N-dimensional SURF descriptor space to estimate the visible body region".
Feulner J, Zhou SK, Angelopoulou E, Seifert S, Cavallaro A, Hornegger J, Comaniciu D.
