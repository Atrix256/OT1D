#pragma once

// y = 1
struct PDFUniform
{
    static const float inline c_xmin = 0.0f;
    static const float inline c_xmax = 1.0f;

    float PDF(float x) const
    {
        if (x < c_xmin || x > c_xmax)
            return 0.0f;
        return 1.0f;
    }

    float ICDF(float x) const
    {
        if (x < c_xmin)
            return 0.0f;

        if (x > c_xmax)
            return 1.0f;

        return x;
    }
};

// y = 2x
struct PDFLinear
{
    static const float inline c_xmin = 0.0f;
    static const float inline c_xmax = 1.0f;

    float PDF(float x) const
    {
        if (x < c_xmin || x > c_xmax)
            return 0.0f;
        return 2.0f * x;
    }

    float ICDF(float x) const
    {
        if (x < c_xmin)
            return 0.0f;

        if (x > c_xmax)
            return 1.0f;

        return std::sqrt(x);
    }
};

// y = 3x^2
struct PDFQuadratic
{
    static const float inline c_xmin = 0.0f;
    static const float inline c_xmax = 1.0f;

    float PDF(float x) const
    {
        if (x < c_xmin || x > c_xmax)
            return 0.0f;
        return 3.0f * x * x;
    }

    float ICDF(float x) const
    {
        if (x < c_xmin)
            return 0.0f;

        if (x > c_xmax)
            return 1.0f;

        return std::powf(x, 1.0f / 3.0f);
    }
};