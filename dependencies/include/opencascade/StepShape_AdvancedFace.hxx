// Created on: 1995-12-01
// Created by: EXPRESS->CDL V0.2 Translator
// Copyright (c) 1995-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#ifndef _StepShape_AdvancedFace_HeaderFile
#define _StepShape_AdvancedFace_HeaderFile

#include "Standard.hxx"
#include "Standard_Type.hxx"

#include "StepShape_FaceSurface.hxx"


class StepShape_AdvancedFace;
DEFINE_STANDARD_HANDLE(StepShape_AdvancedFace, StepShape_FaceSurface)


class StepShape_AdvancedFace : public StepShape_FaceSurface
{

public:

  
  //! Returns a AdvancedFace
  Standard_EXPORT StepShape_AdvancedFace();




  DEFINE_STANDARD_RTTIEXT(StepShape_AdvancedFace,StepShape_FaceSurface)

protected:




private:




};







#endif // _StepShape_AdvancedFace_HeaderFile