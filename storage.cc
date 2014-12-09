
			Mb (X, X) = A (X);
			Mb (Y, X) = A (Y);
			Mb (Z, X) = A (Z);
			Mb (X, Y) = e2 (X);
			Mb (Y, Y) = e2 (Y);
			Mb (Z, Y) = e2 (Z);
			Mb (X, Z) = rW (X);
			Mb (Y, Z) = rW (Y);
			Mb (Z, Z) = rW (Z);
			
			//std::cout<<e1<<e2<<rW<<A<<std::endl;
			//std::cout<<"Mb: "<<Mb<<std::endl;

			float beta = Mb (X, X) * ((Mb(Y, Y) * Mb(Z, Z)) - (Mb (Y, Z) * Mb(Z, Y)))
						- Mb (X, Y) * ((Mb(Z, Z) * Mb(Y, X)) - (Mb (Y, Z) * Mb(Z, X)))
						+ Mb (X, Z) * ((Mb(Y, X) * Mb(Z, Y)) - (Mb (Y, Y) * Mb(Z, X)));
			beta = beta / det;
			//std::cout<<"beta: "<<std::endl;
			
			if (beta < 0.0 || beta > 1.0) 
				break;

			Mg (X, X) = A (X);
			Mg (Y, X) = A (Y);
			Mg (Z, X) = A (Z);
			Mg (X, Y) = e2 (X);
			Mg (Y, Y) = e2 (Y);
			Mg (Z, Y) = e2 (Z);
			Mg (X, Z) = rW (X);
			Mg (Y, Z) = rW (Y);
			Mg (Z, Z) = rW (Z);
			
			float gamma = Mg (X, X) * ((Mg(Y, Y) * Mg(Z, Z)) - (Mg (Y, Z) * Mg(Z, Y)))
						- Mg (X, Y) * ((Mg(Z, Z) * Mg(Y, X)) - (Mg (Y, Z) * Mg(Z, X)))
						+ Mg (X, Z) * ((Mg(Y, X) * Mg(Z, Y)) - (Mg (Y, Y) * Mg(Z, X)));

			gamma = gamma / det;
			if (gamma < 0.0 || beta + gamma > 1.0)
				break;
			
			Mt (X, X) = e1 (X);
			Mt (Y, X) = e1 (Y);
			Mt (Z, X) = e1 (Z);
			Mt (X, Y) = e2 (X);
			Mt (Y, Y) = e2 (Y);
			Mt (Z, Y) = e2 (Z);
			Mt (X, Z) = A (X);
			Mt (Y, Z) = A (Y);
			Mt (Z, Z) = A (Z);
			
			float tstar = Mt (X, X) * ((Mt(Y, Y) * Mt(Z, Z)) - (Mt (Y, Z) * Mt(Z, Y)))
						- Mt (X, Y) * ((Mt(Z, Z) * Mt(Y, X)) - (Mt (Y, Z) * Mt(Z, X)))
						+ Mt (X, Z) * ((Mt(Y, X) * Mt(Z, Y)) - (Mt (Y, Y) * Mt(Z, X)));
			
			tstar = tstar / det;
