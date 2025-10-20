#ifndef GPUFIT_TIME_DOMAIN_LORENTZ_1D_CUH_INCLUDED
#define GPUFIT_TIME_DOMAIN_LORENTZ_1D_CUH_INCLUDED

/* Description of the calculate_time_domain_lorentz_1d function
* ===================================================
* Added by Carlo Bevilacqua on the 10/2025
*
* This function calculates a cosine-modulated exponential decay in the time domain;
* this is equivalent to a Lorentzian in the frequency domain.
* The function isn't expected to be normalized nor of 0-offset : f(tao_ns) = A*exp(-pi*tao_ns*width_GHz)*cos(2*pi*shift_GHz*tao_ns)+offset
*
* This function makes use of the user information data to pass in the
* independent variables (X values) corresponding to the data. The X values
* must be of type REAL.
*
* Note that if no user information is provided, the (X) coordinate of the
* first data value is assumed to be (0.0).  In this case, for a fit size of
* M data points, the (X) coordinates of the data are simply the corresponding
* array index values of the data array, starting from zero.
*
* There are three possibilities regarding the X values:
*
*   No X values provided:
*
*       If no user information is provided, the (X) coordinate of the
*       first data value is assumed to be (0.0).  In this case, for a
*       fit size of M data points, the (X) coordinates of the data are
*       simply the corresponding array index values of the data array,
*       starting from zero.
*
*   X values provided for one fit:
*
*       If the user_info array contains the X values for one fit, then
*       the same X values will be used for all fits.  In this case, the
*       size of the user_info array (in bytes) must equal
*       sizeof(REAL) * n_points.
*
*   Unique X values provided for all fits:
*
*       In this case, the user_info array must contain X values for each
*       fit in the dataset.  In this case, the size of the user_info array
*       (in bytes) must equal sizeof(REAL) * n_points * nfits.
*
* Parameters:
*
* parameters: An input vector of model parameters 
*             p[0]: amplitude 
*             p[1]: shift_GHz
*             p[2]: width_GHz
*             p[3]: offset
*
* n_fits: The number of fits.
*
* n_points: The number of data points per fit.
*
* value: An output vector of model function values.
*
* derivative: An output vector of model function partial derivatives.
*
* point_index: The data point index.
*
* fit_index: The fit index.
*
* chunk_index: The chunk index. Used for indexing of user_info.
*
* user_info: An input vector containing user information.
*
* user_info_size: The size of user_info in bytes.
*
* Calling the calculate_linear1d function
* =======================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_time_domain_lorentz_1d(
    REAL const* parameters,
    int const n_fits,
    int const n_points,
    REAL* value,
    REAL* derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char* user_info,
    std::size_t const user_info_size)
{
    // indices

    REAL* user_info_float = (REAL*)user_info;
    REAL x = 0;
    if (!user_info_float)
    {
        x = point_index;
    }
    else if (user_info_size / sizeof(REAL) == n_points)
    {
        x = user_info_float[point_index];
    }
    else if (user_info_size / sizeof(REAL) > n_points)
    {
        int const chunk_begin = chunk_index * n_fits * n_points;
        int const fit_begin = fit_index * n_points;
        x = user_info_float[chunk_begin + fit_begin + point_index];
    }

	const REAL pi = 3.14159265358979323846264;

    // parameters
    REAL A = parameters[0];
    REAL shift_GHz = parameters[1];
    REAL width_GHz = parameters[2];
    REAL offset = parameters[3];

    // value
    REAL f_exp = exp(-pi * x * width_GHz);
    REAL f_cos = cos(2 * pi * shift_GHz * x);
    REAL f_norm = f_exp * f_cos;
    value[point_index] = A * f_norm + offset;

    // derivatives
    REAL* current_derivatives = derivative + point_index;

    current_derivatives[0 * n_points] = f_norm;                                                     // derivative A
    current_derivatives[1 * n_points] = - A * f_exp * sin(2 * pi * shift_GHz * x) * 2*pi*x;         // derivative shift_GHz
    current_derivatives[2 * n_points] = - A * f_norm * pi*x;                                        // derivative width_GHz
    current_derivatives[3 * n_points] = 1;                                                          // derivative offset
}

#endif