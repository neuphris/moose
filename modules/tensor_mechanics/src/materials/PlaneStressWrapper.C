//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PlaneStressWrapper.h"

registerMooseObject("TensorMechanicsApp", PlaneStressWrapper);

InputParameters
PlaneStressWrapper::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute stress after subtracting inelastic strain increments");
  params.addParam<std::vector<MaterialPropertyName>>("inelastic_strain_names",
                                                     "Names of inelastic strain properties");
  params.addRequiredParam<MaterialName>("stress_calculator", "The material stress calculator");
  params.addRequiredParam<MaterialName>("strain_calculator", "The material strain calculator");
  return params;
}

PlaneStressWrapper::PlaneStressWrapper(
    const InputParameters & parameters)
  : Material(parameters),
  _stress(&getMaterialProperty<RankTwoTensor>("stress")),
  _total_strain(&getMaterialProperty<RankTwoTensor>("total_strain")),
  _total_strain_old(&getMaterialPropertyOld<RankTwoTensor>("total_strain")),
  _out_of_plane_strain(declareProperty<Real>("out_of_plane_strain")),
  _out_of_plane_strain_old(getMaterialPropertyOld<Real>("out_of_plane_strain"))
{
}

void
PlaneStressWrapper::initQpStatefulProperties()
{
   _out_of_plane_strain[_qp] = 0.0;
}

void
PlaneStressWrapper::propagateQpStatefulProperties()
{

  _out_of_plane_strain[_qp] = _out_of_plane_strain_old[_qp];
}


void
PlaneStressWrapper::computeQpStress()
{
}

void
PlaneStressWrapper::computeQpJacobian()
{
}

void
PlaneStressWrapper::initialSetup() // put all of this in the initializer list and delete this function
{
  MaterialName stress_models = getParam<MaterialName>("stress_calculator");
  MaterialBase * stress_model = dynamic_cast<MaterialBase *>(&getMaterialByName(stress_models)); // get materialby stress calculator
  if (stress_model)
    _stress_model = stress_model;
  else
    mooseError("Model " + stress_models + " is not compatible with Materialbase");

  MaterialName strain_models = getParam<MaterialName>("strain_calculator");
  MaterialBase * strain_model = dynamic_cast<MaterialBase *>(&getMaterialByName(strain_models));
  if (strain_model)
    _strain_model = strain_model;
  else
    mooseError("Model " + strain_models + " is not compatible with Materialbase");
}


void
PlaneStressWrapper::computeQpProperties()
{
  _strain_model->resetProperties();
  _strain_model -> _out_of_plane_strain_wrapper = _out_of_plane_strain_old[_qp];
  _strain_model -> _old_plane_strain_wrapper = _out_of_plane_strain_old[_qp];
  _strain_model->computeProperties();
  _stress_model->resetProperties();
  _stress_model->computePropertiesAtQp(_qp);
  Real x1 = _out_of_plane_strain_old[_qp];
  Real out_of_plane_stress = (*_stress)[_qp](2,2);

  _strain_model -> _out_of_plane_strain_wrapper =  _out_of_plane_strain_old[_qp] + 1e-6;
  _strain_model->computeProperties();
  _stress_model->computePropertiesAtQp(_qp);
  Real out_of_plane_stress2 = (*_stress)[_qp](2,2);
  Real x2 = _out_of_plane_strain_old[_qp] + 1e-6;
  Real iter = 0;

  while ((abs(out_of_plane_stress2) > 1e-8))
  {
    Real x_n = x2 - out_of_plane_stress2 * (x2 - x1)/ (out_of_plane_stress2 - out_of_plane_stress);
    _strain_model -> _out_of_plane_strain_wrapper = x_n ;
    _strain_model->computeProperties();
    _stress_model->computePropertiesAtQp(_qp);
    out_of_plane_stress = out_of_plane_stress2;
    out_of_plane_stress2 = (*_stress)[_qp](2,2);
    x1 = x2;
    x2 = x_n;
    iter = iter +1;
  }
  _out_of_plane_strain[_qp] = _strain_model -> _out_of_plane_strain_wrapper;
}
