#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>

#define DETERMINISTIC() false
#include "utils.h"

struct PDFUniform
// y = 1
{
    static const float inline c_xmin = 0.0f;
    static const float inline c_xmax = 1.0f;

    static float PDF(float x)
    {
        if (x < c_xmin || x > c_xmax)
            return 0.0f;
        return 1.0f;
    }

    static float ICDF(float x)
    {
        if (x < c_xmin)
            return 0.0f;

        if (x > c_xmax)
            return 1.0f;

        return x;
    }
};

struct PDFLinear
// y = 2x
{
    static const float inline c_xmin = 0.0f;
    static const float inline c_xmax = 1.0f;

    static float PDF(float x)
    {
        if (x < c_xmin || x > c_xmax)
            return 0.0f;
        return 2.0f * x;
    }

    static float ICDF(float x)
    {
        if (x < c_xmin)
            return 0.0f;

        if (x > c_xmax)
            return 1.0f;

        return std::sqrt(x);
    }
};

struct PDFQuadratic
// y = 3x^2
{
    static const float inline c_xmin = 0.0f;
    static const float inline c_xmax = 1.0f;

    static float PDF(float x)
    {
        if (x < c_xmin || x > c_xmax)
            return 0.0f;
        return 3.0f * x * x;
    }

    static float ICDF(float x)
    {
        if (x < c_xmin)
            return 0.0f;

        if (x > c_xmax)
            return 1.0f;

        return std::powf(x, 1.0f / 3.0f);
    }
};

struct PDFNumeric
// y = (x^3-10x^2+5x+11) / 10.417
{
    static const float inline c_xmin = 0.0f;
    static const float inline c_xmax = 1.0f;

    static const int c_PDFSamples = 1000;
    static const int c_CDFSamples = 100;

    static float PDF(float x)
    {
        if (x < c_xmin || x > c_xmax)
            return 0.0f;
        return (x * x * x - 10.0f * x * x + 5 * x + 11) / 10.417f;
    }

    static float ICDF(float x)
    {
        if (!initialized)
        {
            initialized = true;
            Initialize();
        }

        if (x < c_xmin)
            return 0.0f;

        if (x > c_xmax)
            return 1.0f;

        auto it = std::lower_bound(CDFTable.begin(), CDFTable.end(), x);
        if (it == CDFTable.end())
            return 1.0f;

        int upperIndex = int(it - CDFTable.begin());
        int lowerIndex = std::max(upperIndex - 1, 0);

        float lowerValue = CDFTable[lowerIndex];
        float upperValue = CDFTable[upperIndex];

        float fraction = (x - lowerValue) / (upperValue - lowerValue);

        return (float(lowerIndex) + fraction) / float(c_CDFSamples);
    }

    static void Initialize()
    {
        CDFTable.resize(c_CDFSamples, 0.0f);

        for (int pdfIndex = 0; pdfIndex < c_PDFSamples; ++pdfIndex)
        {
            float x = float(pdfIndex) / float(c_PDFSamples - 1);
            int cdfIndex = Clamp(int(x * float(c_CDFSamples)), 0, c_CDFSamples - 1);
            CDFTable[cdfIndex] += PDF(x) / float(c_PDFSamples);
        }

        for (int cdfIndex = 1; cdfIndex < c_CDFSamples; ++cdfIndex)
            CDFTable[cdfIndex] += CDFTable[cdfIndex - 1];
    }

    static std::vector<float> inline CDFTable;
    static bool inline initialized = false;
};

template <typename PDF1, typename PDF2>
float PWassersteinDistance(float p, int numSamples = 10000000)
{
    // https://www.imagedatascience.com/transport/OTCrashCourse.pdf page 45
    // Integral from 0 to 1 of abs(ICDF1(x) - ICDF2(x))^p
    // Then take the result to ^(1/p)
    pcg32_random_t rng = GetRNG();
    double ret = 0.0;
    for (int i = 0; i < numSamples; ++i)
    {
        float x = RandomFloat01(rng);
        float icdf1 = PDF1::ICDF(x);
        float icdf2 = PDF2::ICDF(x);
        double y = std::pow(std::abs((double)icdf1 - (double)icdf2), p);
        ret = Lerp(ret, y, 1.0 / double(i + 1));
    }
    return (float)std::pow(ret, 1.0f / p);
}

int main(int argc, char** argv)
{
    //float d1 = PWassersteinDistance<PDFUniform, PDFUniform>(1.0f);
    //printf("%f\n", d1);

    printf("Uniform To Linear = %f\n", PWassersteinDistance<PDFUniform, PDFLinear>(2.0f));
    printf("Uniform To Quadratic = %f\n", PWassersteinDistance<PDFUniform, PDFQuadratic>(2.0f));
    printf("Linear To Quadratic = %f\n", PWassersteinDistance<PDFLinear, PDFQuadratic>(2.0f));

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