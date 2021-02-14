# Physical-Modeling-Piano-Synthesis

The goal of this project is the development of a virtual piano via matlab.
The approach used is that of physical modeling, that is, the study of interactions that generate a certain sound in the real world.
The sound is obtained by a simulation consisting of three steps, aimed at replicating the mechanics of a real piano:
-hammer-string interaction
-model of the string
-amplification of the oscillation via soundboard

The first phase was implemented using the Finite Differences technique, i.e. discretizing the equations that shape this interaction proposed in askenfelt research, 1990.
The second part of the model implements wave propagation in a piano string using the waveguides technique.  The differential equation describing this phenomenon takes into account dispersion and rigidity of the string (Bensa, 2003).
The propagation of vibration propagating in the soundboard has been approximated here with a convolution between the signal produced at the previous steps and the response to the pulse of a properly recorded grand piano. The latter does not reflect the physical behaviour of the instrument but allows to simplify the model and obtain a piano sound very similar to the real one, but at the expense of a greater demand for resources.

In order to make the sound more realistic, reverberation has also been implemented, allowing the user to choose between three different reverberating environments to shape the sound produced by the previous model.
