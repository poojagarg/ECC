#!/bin/bash
#Run the following commands in the terminal with current directory as "examples"
  
chmod +x *.sh
  
./FiniteFieldArithmetic.sh
./FiniteFieldArithmetic < FiniteFieldArithmeticI
  
./ECarithmetic.sh
./ECarithmetic < ECarithmeticI
./ECarithmetic < ECarithmetic_2_1_I
./ECarithmetic < ECarithmetic_2_2_I

./pohlig.sh
./pohlig < pohligI
./pohlig < pohlig_2_1_I
./pohlig < pohlig_2_2_I

./ELGAMAL_DSA.sh
./ELGAMAL_DSA < ELGAMAL_DSAI
  
./ECDSA.sh
./ECDSA < ECDSAI
  
./ELGAMAL_Encryption.sh
./ELGAMAL_Encryption < ELGAMAL_EncryptionI

#To check examples program for Elliptic Curve of type 1 and type 2, we use input file that ends with _2_1_I and _2_2_I for type 1 and type 2 respectively.

 