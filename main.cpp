#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>

#define DETERMINISTIC() false
#include "utils.h"

#include "analytic.h"
#include "numeric.h"

// y = (x^3-10x^2+5x+11) / 10.417
// return (x * x * x - 10.0f * x * x + 5 * x + 11) / 10.417f;

template <typename PDF1, typename PDF2>
float PWassersteinDistance(float p, const PDF1& pdf1, const PDF2& pdf2, int numSamples = 10000000)
{
    // https://www.imagedatascience.com/transport/OTCrashCourse.pdf page 45
    // Integral from 0 to 1 of abs(ICDF1(x) - ICDF2(x))^p
    // Then take the result to ^(1/p)
    pcg32_random_t rng = GetRNG();
    double ret = 0.0;
    for (int i = 0; i < numSamples; ++i)
    {
        float x = RandomFloat01(rng);
        float icdf1 = pdf1.ICDF(x);
        float icdf2 = pdf2.ICDF(x);
        double y = std::pow(std::abs((double)icdf1 - (double)icdf2), p);
        ret = Lerp(ret, y, 1.0 / double(i + 1));
    }
    return (float)std::pow(ret, 1.0f / p);
}

int main(int argc, char** argv)
{
    //float d1 = PWassersteinDistance<PDFUniform, PDFUniform>(1.0f);
    //printf("%f\n", d1);

    printf("Uniform To Linear = %f\n", PWassersteinDistance(2.0f, PDFUniform(), PDFLinear()));
    printf("Uniform To Quadratic = %f\n", PWassersteinDistance(2.0f, PDFUniform(), PDFQuadratic()));
    printf("Linear To Quadratic = %f\n", PWassersteinDistance<>(2.0f, PDFLinear(), PDFQuadratic()));

    //float d1 = PWassersteinDistance<PDFUniform, PDFLinear>(1.0f);
    //float d2 = PWassersteinDistance<PDFUniform, PDFLinear>(2.0f);
    //printf("%f\n%f\n", d1, d2);

    return 0;
}

/*
TODO:
- interpolate PDFs. show by drawing random numbers from it.
 ? i think this is by lerping a CDF histogram and interpolating, then drawing from the CDF.

? how to do interpolation between PDFs?
 * both numerical and analytical.

* do some wasserstein distance analytically (in blog post, but also in general)
 * woohoo! integration of polynomials. Can integrate each term individually. very easy.
 * probably need this for analytical interpolation between PDFs.

NOTES:
- links from email
- link to this for inverted CDF talk: https://blog.demofox.org/2017/08/05/generating-random-numbers-from-a-specific-distribution-by-inverting-the-cdf/
- this works with tabulated CDFs as well, doesn't have to be analytical.
*/