#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>

template <typename T>
T Clamp(T x, T themin, T themax)
{
    if (x <= themin)
        return themin;
    if (x >= themax)
        return themax;
    return x;
}

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

int main(int argc, char** argv)
{
    float a = PDFNumeric::ICDF(0.5f);
    float b = PDFNumeric::PDF(a);

    return 0;
}

/*
TODO:
- P-wasserstein where people can be 1 or 2 or whatever.
- interpolate PDFs. show by drawing random numbers from it.
- measure difference between two pdfs?
- analytical and tabular or just analytical?


NOTES:
- links from email
- link to this for inverted CDF talk: https://blog.demofox.org/2017/08/05/generating-random-numbers-from-a-specific-distribution-by-inverting-the-cdf/
- this works with tabulated CDFs as well, doesn't have to be analytical.
*/