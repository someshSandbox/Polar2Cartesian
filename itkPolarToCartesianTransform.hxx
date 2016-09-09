/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkPolarToCartesianTransform_hxx
#define itkPolarToCartesianTransform_hxx

#include "itkPolarToCartesianTransform.h"

namespace itk
{
// Constructor with default arguments
template<typename TParametersValueType, unsigned int NDimensions>
PolarToCartesianTransform<TParametersValueType, NDimensions>
::PolarToCartesianTransform()
// add this construction call when deriving from itk::Transform
// :Superclass(ParametersDimension)
{
  m_ForwardPolarToCartesian = true;
}

// Destructor
template<typename TParametersValueType, unsigned int NDimensions>
PolarToCartesianTransform<TParametersValueType, NDimensions>::
~PolarToCartesianTransform()
{
}

// Print self
template<typename TParametersValueType, unsigned int NDimensions>
void
PolarToCartesianTransform<TParametersValueType, NDimensions>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "x = r * std::cos(theta) " << std::endl;
  os << indent << "y = r * std::sin(theta) " << std::endl;
  os << indent << "r = sqrt(x*x + y*y)" << std::endl;
  os << indent << "theta = atan(y/x)" << std::endl;
  os << indent << "m_ForwardPolarToCartesian = " << ( m_ForwardPolarToCartesian ? "True" : "False" );
  os << indent << std::endl;
}

template<typename TParametersValueType, unsigned int NDimensions>
typename PolarToCartesianTransform<TParametersValueType, NDimensions>
::OutputPointType
PolarToCartesianTransform<TParametersValueType, NDimensions>::TransformPoint(const InputPointType & point) const
{
  OutputPointType result;

  if ( m_ForwardPolarToCartesian )
  {
    result = TransformPolarToCartesian(point);
  }
  else
  {
    result = TransformCartesianToPolar(point);
  }
  return result;
}

/** Transform a point, from cylindrical to cartesian */
template<typename TParametersValueType, unsigned int NDimensions>
typename PolarToCartesianTransform<TParametersValueType, NDimensions>
::OutputPointType
PolarToCartesianTransform<TParametersValueType,
                                NDimensions >
                                ::TransformPolarToCartesian(const InputPointType & point) const
{

  // point[0] = R
  // point[1] = Theta
 
  ScalarType degreeToRadian = vnl_math::pi / 180.0;
  ScalarType r              = point[0];
  ScalarType theta          = degreeToRadian * point[1];

 OutputPointType result = point; // Converted point
   // X
  result[0] = r * std::cos(theta);
  // Y
  result[1] = r * std::sin(theta);
 
  return result;
}

template<typename TParametersValueType, unsigned int NDimensions>
typename PolarToCartesianTransform<TParametersValueType, NDimensions>
::OutputPointType
PolarToCartesianTransform<TParametersValueType, NDimensions>::TransformCartesianToPolar(
  const OutputPointType & point) const
{
  // point[0] = X
  // point[1] = Y
 OutputPointType result = point;       // Converted point
  ScalarType radianToDegrees = 180.0 / vnl_math::pi ;
  // R
  result[0] = std::sqrt(point[0] * point[0]  + 
                        point[1] * point[1]);
  // Theta
  result[1] = std::atan2(point[1], point[0]) * radianToDegrees;  // Theta
  if(result[1] < 0)
    result[1] += 360.0;
 
  return result;
}



template<typename TParametersValueType, unsigned int NDimensions>
void
PolarToCartesianTransform<TParametersValueType, NDimensions>::SetForwardPolarToCartesian()
{
  m_ForwardPolarToCartesian = true;
}

template<typename TParametersValueType, unsigned int NDimensions>
void
PolarToCartesianTransform<TParametersValueType, NDimensions>::SetForwardCartesianToPolar()
{
  m_ForwardPolarToCartesian = false;
}
} //namespace
#endif
