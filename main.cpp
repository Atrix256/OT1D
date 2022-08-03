#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>

#define DETERMINISTIC() false
#include "utils.h"

#include "analytic.h"
#include "numeric.h"

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

template <typename PDF1, typename PDF2>
void InterpolatePDFs_PDF(const char* fileName, const PDF1& pdf1, const PDF2& pdf2, int numSteps = 5, int numValues = 100)
{
    printf("%s...\n", fileName);

    // Make the interpolated PDFs
    std::vector<std::vector<float>> PDFs(numSteps);
    for (int step = 0; step < numSteps; ++step)
    {
        // Make the PDF
        float t = float(step) / float(numSteps - 1);
        std::vector<float>& PDF = PDFs[step];
        PDF.resize(numValues, 0.0f);
        for (int i = 0; i < numValues; ++i)
        {
            float x = float(i) / float(numValues - 1);
            float y1 = pdf1.PDF(x);
            float y2 = pdf2.PDF(x);
            PDF[i] = Lerp(y1, y2, t);
        }

        // normalize PDF
        float total = 0.0f;
        for (float f : PDF)
            total += f;
        for (float& f : PDF)
            f /= total;
    }

    // Write it to a file
    FILE* file = nullptr;
    fopen_s(&file, fileName, "wt");

    for (int column = 0; column < numSteps; ++column)
        fprintf(file, "\"t=%i%%\",", int(100.0f * float(column) / float(numSteps - 1)));
    fprintf(file, "\n");

    for (int row = 0; row < numValues; ++row)
    {
        for (int column = 0; column < numSteps; ++column)
            fprintf(file, "\"%f\",", PDFs[column][row]);
        fprintf(file, "\n");
    }

    for (int column = 0; column < numSteps; ++column)
    {
        float total = 0.0f;
        for (float f : PDFs[column])
            total += f;

        printf("Column %i total = %0.2f\n", column, total);
    }

    fclose(file);
    printf("\n");
}

// TODO: tune numValuesICDF

template <typename PDF1, typename PDF2>
void InterpolatePDFs_ICDF(const char* fileName, const PDF1& pdf1, const PDF2& pdf2, int numSteps = 5, int numValuesICDF = 1000000, int numValuesPDF = 100)
{
    printf("%s...\n", fileName);

    // Make the interpolated PDFs
    std::vector<std::vector<float>> PDFs(numSteps);
    for (int step = 0; step < numSteps; ++step)
    {
        // make the ICDF
        float t = float(step) / float(numSteps - 1);
        std::vector<float> ICDF;
        ICDF.resize(numValuesICDF, 0.0f);
        for (int i = 0; i < numValuesICDF; ++i)
        {
            float x = float(i) / float(numValuesICDF - 1);
            float y1 = pdf1.ICDF(x);
            float y2 = pdf2.ICDF(x);
            ICDF[i] = Lerp(y1, y2, t);
        }
        ICDF[numValuesICDF - 1] = 1.0f;

        // TODO: CDF made PDF looks 1 or half a pixel shifted to the left.

        // make the CDF by inverting the ICDF
        std::vector<float> CDF;
        CDF.resize(numValuesPDF + 1, 0.0f);
        for (int i = 0; i <= numValuesPDF; ++i)
        {
            float x = float(i) / float(numValuesPDF);

            auto it = std::lower_bound(ICDF.begin(), ICDF.end(), x);
            if (it == ICDF.end())
            {
                printf("Could not find value %f in ICDF table! (index %i/%i)\n", x, i, numValuesPDF);
            }
            else
            {
                int upperIndex = int(it - ICDF.begin());
                int lowerIndex = std::max(upperIndex - 1, 0);

                if (upperIndex == lowerIndex)
                {
                    CDF[i] = float(lowerIndex) / float(numValuesPDF);
                }
                else
                {
                    float lowerValue = ICDF[lowerIndex];
                    float upperValue = ICDF[upperIndex];

                    float fraction = (x - lowerValue) / (upperValue - lowerValue);

                    CDF[i] = (float(lowerIndex) + fraction) / float(numValuesPDF);
                }
            }
        }

        // normalize the CDF
        for (float& f : CDF)
            f /= CDF[numValuesPDF];

        // make the PDF from the CDF
        std::vector<float>& PDF = PDFs[step];
        PDF.resize(numValuesPDF, 0.0f);
        for (int i = 0; i < numValuesPDF; ++i)
            PDF[i] = CDF[i + 1] - CDF[i];

        // normalize the PDF
        float total = 0.0f;
        for (float f : PDF)
            total += f;
        for (float& f : PDF)
            f /= total;
    }

    // Write it to a file
    FILE* file = nullptr;
    fopen_s(&file, fileName, "wt");

    for (int column = 0; column < numSteps; ++column)
        fprintf(file, "\"t=%i%%\",", int(100.0f * float(column) / float(numSteps - 1)));
    fprintf(file, "\n");

    for (int row = 0; row < numValuesPDF; ++row)
    {
        for (int column = 0; column < numSteps; ++column)
            fprintf(file, "\"%f\",", PDFs[column][row]);
        fprintf(file, "\n");
    }

    fclose(file);

    printf("\n");
}

int main(int argc, char** argv)
{
#if 0
    PDFNumeric pdftTableUniform([](float x) { return 1.0f; });
    PDFNumeric pdftTableLinear([](float x) { return 2.0f * x; });
    PDFNumeric pdftTableQuadratic([](float x) { return 3.0f * x * x; });

    printf("(analytical p=2) Uniform To Linear = %f\n", PWassersteinDistance(2.0f, PDFUniform(), PDFLinear()));
    printf("(analytical p=2) Uniform To Quadratic = %f\n", PWassersteinDistance(2.0f, PDFUniform(), PDFQuadratic()));
    printf("(analytical p=2) Linear To Quadratic = %f\n\n", PWassersteinDistance<>(2.0f, PDFLinear(), PDFQuadratic()));

    printf("(table p=2) Uniform To Linear = %f\n", PWassersteinDistance(2.0f, pdftTableUniform, pdftTableLinear));
    printf("(table p=2) Uniform To Quadratic = %f\n", PWassersteinDistance(2.0f, pdftTableUniform, pdftTableQuadratic));
    printf("(table p=2) Linear To Quadratic = %f\n\n", PWassersteinDistance(2.0f, pdftTableLinear, pdftTableQuadratic));

    printf("(table p=1) Uniform To Linear = %f\n", PWassersteinDistance(1.0f, pdftTableUniform, pdftTableLinear));
    printf("(table p=1) Uniform To Quadratic = %f\n", PWassersteinDistance(1.0f, pdftTableUniform, pdftTableQuadratic));
    printf("(table p=1) Linear To Quadratic = %f\n\n", PWassersteinDistance(1.0f, pdftTableLinear, pdftTableQuadratic));

    printf("(table p=3) Uniform To Linear = %f\n", PWassersteinDistance(3.0f, pdftTableUniform, pdftTableLinear));
    printf("(table p=3) Uniform To Quadratic = %f\n", PWassersteinDistance(3.0f, pdftTableUniform, pdftTableQuadratic));
    printf("(table p=3) Linear To Quadratic = %f\n\n", PWassersteinDistance(3.0f, pdftTableLinear, pdftTableQuadratic));
#endif

    PDFNumeric pdftTableGauss1([](float x) { x -= 0.2f; return exp(-x * x / (2.0f * 0.1f * 0.1f)); });
    PDFNumeric pdftTableGauss2([](float x) { x -= 0.6f; return exp(-x * x / (2.0f * 0.15f * 0.15f)); });
    InterpolatePDFs_PDF("_Gauss2Gauss_PDF.csv", pdftTableGauss1, pdftTableGauss2);
    InterpolatePDFs_ICDF("_Gauss2Gauss_CDF.csv", pdftTableGauss1, pdftTableGauss2);

    // TODO: do other p values instead of only 1.
    // TODO: do 3 way interpolation?

    return 0;
}

/*
TODO:
- interpolate PDFs. show by drawing random numbers from it.
 ? i think this is by lerping a CDF histogram and interpolating, then drawing from the CDF.
 * try it in other P's.
 * NOTE: i guess the difference is interpolating CDFs (or ICDFs?) instead of PDFs?
 * NOTE: link to bezier curves and bezier triangles as other interpolation forms

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