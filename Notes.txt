Envelope Space Charge Files:


SC_UprightEnv.h	      Header file for upright (non-tilted, KV) envelope
		      solver.
	_KickEnvelope   - Applies space charge terms to envelope, which
			  is stored in first macroparticle.
	_BoundaryCoeffs - Calculates coefficients for a conducting wall
			  boundary if one is specified.
	_ApplyForce	- Applies space charge terms to macroparticles
			  calculated according to the envelope parameters.
			  The macroparticles are stored starting in the
			  second location.
	_getu		- Calculates auxilliary parameter u in space
			  charge formulation when particle is outside the
			  envelope.
	_gettheta	- Calculates auxilliary parameter theta in space
			  charge formulation when particle is outside the
			  envelope.

SC_UprightEnv.cc      c++ code for upright (non-tilted, KV) envelope solver.


SC_TiltEnv.h	      Header file for tilted (rotating beam) envelope solver.
	_KickEnvelope   - Applies space charge terms to envelope, which
			  is stored in first and second macroparticles.
	_BoundaryCoeffs - Calculates coefficients for a conducting wall
			  boundary if one is specified.
	_ApplyForce	- Applies space charge terms to macroparticles
			  calculated according to the envelope parameters.
			  The macroparticles are stored starting in the
			  third location.
	_getu		- Calculates auxilliary parameter u in space
			  charge formulation when particle is outside the
			  envelope.
	_gettheta	- Calculates auxilliary parameter theta in space
			  charge formulation when particle is outside the
			  envelope.
	_getEllipsDat	- Calculates major and minor axes and tilt of
			  ellipse given the envelope parameters. It will
			  work for any phase plane ellipse if written in
			  the Danilov parameterization.


SC_TiltEnv.cc	      c++ code for tilted (rotating beam) envelope solver.


TSpaceCharge.cc	      c++ code from the old ORBIT Code illustrating the use
		      of the upright and tilted solvers in the ORBIT Code.
