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
#ifndef itkPolarToCartesianTransform_h
#define itkPolarToCartesianTransform_h

#include "itkAffineTransform.h"
#include "vnl/vnl_math.h"

namespace itk
{
/** \class PolarToCartesianTransform
 * \brief Transforms from an cylindrical coordinate system to
 * a Cartesian coordinate system, or vice versa.
 *
 * The three coordinate axis are r, theta, and z.
 *
 * The equations form performing the conversion from cylindrical
 * coordinates to cartesian coordinates are as follows:
 * x = r * std::cos(theta)
 * y = r * std::sin(theta)
 * z = z
 *
 * The reversed transforms are:
 * r = sqrt(x*x + y*y)
 * theta = atan(y/x)
 * z = z
 *
 * In this class, we can also set what a "forward" transform means.  If we call
 * SetForwardPolarToCartesian(), a forward transform will return
 * cartesian coordinates when passed cylindrical coordinates.  Calling
 * SetForwardCartesianToPolar() will cause the forward transform
 * to return cylindricalcoordinates from cartesian coordinates.
 *
 * Setting the FirstSampleDistance to a non-zero value means that a r value
 * of 12 is actually (12 + FirstSampleDistance) distance from the origin.
 *
 * There are two template parameters for this class:
 *
 * ScalarT       The type to be used for scalar numeric values.  Either
 *               float or double.
 *
 * NDimensions   The number of dimensions of the vector space (must be >=3).
 *
 * \todo Is there any real value in allowing the user to template
 * over the scalar type?  Perhaps it should always be double, unless
 * there's a compatibility problem with the Point class.
 *
 * \todo  Derive this class from a yet undefined TransformBase class.
 *        Currently, this class derives from AffineTransform, although
 *        it is not an affine transform.
 * \ingroup ITKTransform
 *
 * \wiki
 * \wikiexample{Utilities/PolarToCartesianTransform,Cartesian to Polar and vice-versa}
 * \endwiki
 */
template<typename TParametersValueType = double,
         unsigned int NDimensions = 3>
class PolarToCartesianTransform:
  public AffineTransform<TParametersValueType, NDimensions>
{
public:
  /** Standard class typedefs.   */
  typedef PolarToCartesianTransform               Self;
  typedef AffineTransform<TParametersValueType, NDimensions> Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;

  /** Dimension of the domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
  itkStaticConstMacro( ParametersDimension, unsigned int,
                       NDimensions * ( NDimensions + 1 ) );

  /** Run-time type information (and related methods).   */
  itkTypeMacro(PolarToCartesianTransform, AffineTransform);

  /** New macro for creation of through a Smart Pointer.   */
  itkNewMacro(Self);

  /** Parameters type.   */
  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::FixedParametersType FixedParametersType;

  /** Jacobian type.   */
  typedef typename Superclass::JacobianType JacobianType;

  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard coordinate point type for this class */
  typedef  typename Superclass::InputPointType  InputPointType;
  typedef  typename Superclass::OutputPointType OutputPointType;

  /** Standard matrix type for this class.   */
  typedef Matrix< TParametersValueType, itkGetStaticConstMacro(SpaceDimension),
          itkGetStaticConstMacro(SpaceDimension) > MatrixType;


  /** Transform from azimuth-elevation to cartesian. */
  OutputPointType     TransformPoint(const InputPointType  & point) const ITK_OVERRIDE;

  /** Back transform from cartesian to azimuth-elevation.  */
  inline InputPointType  BackTransform(const OutputPointType  & point) const
  {
    InputPointType result;

    if ( m_ForwardPolarToCartesian )
    {
      result = static_cast< InputPointType >( TransformCartesianToPolar(point) );
    }
    else
    {
      result = static_cast< InputPointType >( TransformPolarToCartesian(point) );
    }
    return result;
  }

  inline InputPointType  BackTransformPoint(const OutputPointType  & point) const
  {
    return BackTransform(point);
  }

  /** Defines that the forward transform goes from azimuth,elevation to
   *  cartesian. */
  void SetForwardPolarToCartesian();

  /** Defines that the forward transform goes from cartesian to azimuth,
   *  elevation.  */
  void SetForwardCartesianToPolar();

  /** Perform conversion from Azimuth Elevation coordinates to Cartesian
   *  Coordinates. */
  OutputPointType TransformPolarToCartesian(const InputPointType & point) const;

  /** Perform conversion from Cartesian Coordinates to Azimuth Elevation
   *  coordinates.  */
  OutputPointType TransformCartesianToPolar(const OutputPointType & point) const;


protected:
  /** Create an PolarToCartesianTransform object. */
  PolarToCartesianTransform();

  /** Destroy an PolarToCartesianTransform object. */
  virtual ~PolarToCartesianTransform();

  /** Print contents of an PolarTransform. */
  void PrintSelf(std::ostream & s, Indent indent) const ITK_OVERRIDE;

private:
  PolarToCartesianTransform(const Self &); // purposely not
  // implemented
  void operator=(const Self &);                       //purposely not

  // implemented

  bool   m_ForwardPolarToCartesian;

}; //class PolarToCartesianTransform
}  // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPolarToCartesianTransform.hxx"
#endif

#endif /* itkPolarToCartesianTransform_h */
