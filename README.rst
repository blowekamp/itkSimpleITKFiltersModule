General
=======

.. image:: https://ci.appveyor.com/api/projects/status/je74aykgrftxxkdd/branch/master?svg=true     
    :target: https://ci.appveyor.com/project/blowekamp/itksimpleitkfiltersmodule-muow0/branch/master
    
.. image:: https://travis-ci.org/SimpleITK/itkSimpleITKFiltersModule.svg?branch=master
    :target: https://travis-ci.org/SimpleITK/itkSimpleITKFiltersModule
     
This is a module for `ITK <http://itk.org>`_: The Insight Toolkit for Segmentation and
Registration. It is designed to work with the ITKv4
modular system and  to be places it ITK/Modules/External or uses as a
Remote module.

This SimpleITKFilters module is a collection of ITK Filters used by SimpleITK.

Getting Started
===============

The official ITKv4 Wiki documentation on adding an external module is here:

* http://www.itk.org/Wiki/ITK_Release_4/Modularization/Add_an_external_module_(external_module)
* http://www.itk.org/Wiki/ITK/Policy_and_Procedures_for_Adding_Remote_Modules


External Module
---------------

The following is a breif list of instructions to get a external module
started in a repository::

  mkdir ITK/Modules/External/itkMyModule
  cd ITK/Modules/External/itkMyModule
  git init
  git fetch git@github.com:blowekamp/itkExternalTemplate.git
  git merge FETCH_HEAD

Remote Module
-------------

After a module has been created as a git repository it can be included
as a remote module, which enables automatic fetching. Add a file in
"ITK/Modules/Remote" called "YourModule.remote.cmake", for this module
it would be "ExternalExample.remote.cmake" with the followlowing contents::

  itk_fetch_module(ExternalTemplate
    "A template for a module."
    GIT_REPOSITORY http://github.com/blowekamp/itkExternalTemplate.git
    GIT_TAG master
    )

License
=======

This software is distributed under the Apache License. Please see
LICENSE for details.


Authors
=======

* Bradley Lowekamp
