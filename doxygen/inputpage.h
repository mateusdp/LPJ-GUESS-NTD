/** \page input_page Input modules - providing forcing data to the model

LPJ-GUESS gets its forcing data from an input module, a subclass of the InputModule class.
The model is distributed with the input modules CRUInput (for the CRU dataset) and DemoInput
(for demonstration purposes, reads in a toy dataset).

The user of the model will often need to either modify one of the supplied input modules or
create a new input module. In some cases it may also be possible to inherit from an existing
input module and adjust their behaviour by overriding some of its member functions.

- \ref InputModule
- \ref CRUInput
- \ref DemoInput
- \ref GetclimInput
- \ref REGISTER_INPUT_MODULE - A macro for registering a new input module
*/
