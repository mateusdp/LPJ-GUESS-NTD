/** \page input_page Input modules - providing forcing data to the model

Climate forcing input:

LPJ-GUESS gets its climate-forcing data from an input module, a subclass of the InputModule class.
The model is distributed with the input modules CRUInput (for the CRU dataset), DemoInput
(for demonstration purposes, reads in a toy dataset), CFInput (a simple example implementation for 
reading NetCDF files), and GetclimInput

The user of the model will often need to either modify one of the supplied input modules or
create a new input module. In some cases it may also be possible to inherit from an existing
input module and adjust their behaviour by overriding some of its member functions.

Base class
- \ref InputModule

Climate forcing.
- \ref CRUInput
- \ref DemoInput
- \ref CFInput
- \ref GetclimInput

A macro for registering a new input module
- \ref REGISTER_INPUT_MODULE


Other input classes:

Input classes that do not inherit from InputModule.

- \ref LandcoverInput
- \ref ManagementInput
- \ref SoilInput

*/
