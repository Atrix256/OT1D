#pragma once

#include <vector>
#include <functional>

struct PDFNumeric
{
    static const float inline c_xmin = 0.0f;
    static const float inline c_xmax = 1.0f;

    static const int c_PDFSamples = 10000;
    static const int c_CDFSamples = 100;

    typedef std::function<float(float)> PDFFn;

    PDFNumeric(PDFFn pdf)
    {
        m_PDF = pdf;

        // Make a discretized PDF table
        m_CDFTable.resize(c_CDFSamples, 0.0f);
        for (int pdfIndex = 0; pdfIndex < c_PDFSamples; ++pdfIndex)
        {
            float x = float(pdfIndex) / float(c_PDFSamples - 1);
            int cdfIndex = Clamp(int(x * float(c_CDFSamples)), 0, c_CDFSamples - 1);
            m_CDFTable[cdfIndex] += PDF(x);
        }

        // normalize PDF to sum to 1.0
        float total = 0.0f;
        for (float f : m_CDFTable)
            total += f;
        for (float& f : m_CDFTable)
            f /= total;

        // Make CDF table
        for (int cdfIndex = 1; cdfIndex < c_CDFSamples; ++cdfIndex)
            m_CDFTable[cdfIndex] += m_CDFTable[cdfIndex - 1];

        // normalize CDF so last value is 1.0
        for (float& f : m_CDFTable)
            f /= m_CDFTable[c_CDFSamples - 1];
    }

    float PDF(float x) const
    {
        if (x < c_xmin || x > c_xmax)
            return 0.0f;
        return m_PDF(x);
    }

    float CDF(float x) const
    {
        if (x < c_xmin)
            return 0.0f;

        if (x > c_xmax)
            return 1.0f;

        float index = Clamp(x * float(c_CDFSamples), 0.0f, float(c_CDFSamples - 1));

        int index1 = int(index);
        int index2 = Clamp(index1 + 1, 0, c_CDFSamples - 1);
        float fract = index - floor(index);
        return Lerp(m_CDFTable[index1], m_CDFTable[index2], fract);
    }

    float ICDF(float x) const
    {
        if (x < c_xmin)
            return 0.0f;

        if (x > c_xmax)
            return 1.0f;

        auto it = std::lower_bound(m_CDFTable.begin(), m_CDFTable.end(), x);
        if (it == m_CDFTable.end())
            return 1.0f;

        int upperIndex = int(it - m_CDFTable.begin());
        int lowerIndex = std::max(upperIndex - 1, 0);

        if (lowerIndex == upperIndex)
            return float(lowerIndex) / float(c_CDFSamples);

        float lowerValue = m_CDFTable[lowerIndex];
        float upperValue = m_CDFTable[upperIndex];

        float fraction = (x - lowerValue) / (upperValue - lowerValue);

        return (float(lowerIndex) + fraction) / float(c_CDFSamples);
    }

    PDFFn m_PDF;
    std::vector<float> m_CDFTable;
};