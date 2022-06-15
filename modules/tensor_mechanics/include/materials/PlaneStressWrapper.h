//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

class PlaneStressWrapper : public Material
{
public:
  static InputParameters validParams();

  PlaneStressWrapper(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() ;
  virtual void computeQpStress();
  virtual void computeQpJacobian();
  virtual void initialSetup() override;
  virtual void computeQpProperties() override;

  /// Material property for current stress
  const MaterialProperty<RankTwoTensor> * _stress;

  // /// Material property for current strain
   const MaterialProperty<RankTwoTensor> * _total_strain;

  /// Material property for old strain
  const MaterialProperty<RankTwoTensor> * _total_strain_old;

  MaterialBase * _stress_model;
  MaterialBase * _strain_model;

  MaterialProperty<Real> & _out_of_plane_strain;
  const MaterialProperty<Real> & _out_of_plane_strain_old;

};
