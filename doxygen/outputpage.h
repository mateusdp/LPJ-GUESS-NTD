/** \page output_page Output modules - generating output files with model results

LPJ-GUESS is distributed with an output module for many common outputs (CommonOutput), 
but the user may need to supply code for additional outputs. This can be done either 
by modifying the standard output module, or by creating a new output module.

All output modules inherit from the base class OutputModule.

Most output related code is gathered in the namespace GuessOutput.

- \ref GuessOutput::CommonOutput
- \ref GuessOutput::MiscOutput
- \ref GuessOutput::OutputModule
- \ref REGISTER_OUTPUT_MODULE - A macro for registering a new output module

The output modules typically generate text output. They should however not create these
text files directly, but rather use convenience classes. All output should go through a
so called output channel. The output channel will take care of generating the correct file
format, and for instance align columns etc.

The output channel related classes are defined in outputchannel.h.
*/
