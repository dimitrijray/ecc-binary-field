#Importing modules
import math
import binop
import random
from primes import Primes

# ====== Global variables and functions ======
# -- Constants

# Point at infinity
O = ['0','0']

# -- Functions

# Fast exponentiation
def fastexp(n,exp):
    result = 1
    multi = n
    while exp>0:
        if(exp%2==1):
            result = result*multi
        exp=exp>>1
        multi = multi**2
    return result

# Create an elliptic curve object.
class Curve(object):
    # The __init__ module generates the curve and its parameters.
    #-------------------------------------------------------------
    def __init__(self,deg,manual = False,a = '0', b = '0'):
        self.deg = deg
        #Calculate the appropriate irreducible polynomial for this curve.
        self.mod = binop.genIrred(deg)
        #Check whether the manual mode is set or not.
        ## If the manual mode is set, take the 'a' and 'b' input as
        ## the curve parameters.  Else, randomize it.
        if (manual):
            self.a = a
            self.b = b
        else:
            temp = random.randrange(1,(2**deg))
            self.a = bin(temp)[2:]
            temp = random.randrange(1,(2**deg))
            self.b = bin(temp)[2:]

    # Point generation for the ElGamal
    #-------------------------------------------------------------
    def generatePoint(self):
        isPoint = False
        # For debugging purposes: calculating no. of iterations.
        nIterate = 1
        print("Generating points...")
        while(not isPoint):
            # Select a random X - coordinate.
            temp = random.randrange(1,fastexp(2,self.deg))
            X = bin(temp)[2:]
            # Calculate the RHS, which depends solely on X.
            RHS = binop.add(binop.add(binop.powMod(X,3,self.mod),binop.mod(binop.multi(self.a,binop.square(X)),self.mod)),self.b)
            # Divide RHS with X^(-2)
            B = binop.mod(binop.square(X),self.mod)
            Sol = binop.divide(RHS,B,self.mod)
            # Calculate the trace of Sol.
            Tr = '0'
            adder = Sol
            for e in range(self.deg):
                trexp = fastexp(2,e)
                Tr = binop.add(Tr,adder)
                adder = binop.mod(binop.square(adder),self.mod)
            # If the Trace is zero, then the selected X-coordinate has a solution in Y.  Find the solution
            # using Kugurakov's Formula
            # Else, go find another X-coordinate.
            if (Tr == '0'):
                #Find a "Y-coordinate" with a trace of 1.
                isSolved = False
                while(not isSolved):
                    temp = random.randrange(1,fastexp(2,self.deg))
                    C  = bin(temp)[2:]
                    TrC = '0'
                    adderC = C
                    for e in range(self.deg):
                        trexpc = fastexp(2,e)
                        TrC = binop.add(TrC,adderC)
                        adderC = binop.mod(binop.square(adderC),self.mod)
                    if (TrC == '1'):
                        #If an element with trace 1 is found, find the solution!
                        isSolved = True
                        Y = '0'
                        delta = binop.mod(binop.square(Sol),self.mod)
                        for j in range(1,self.deg):
                            adderY = '0'
                            u = C
                            for k in range(j):
                                adderY = binop.add(adderY,u)
                                u = binop.mod(binop.square(u),self.mod)
                            adderY = binop.mod(binop.multi(delta,adderY),self.mod)
                            Y = binop.add(Y,adderY)
                            delta = binop.mod(binop.square(delta),self.mod)
                #Done with the solution. Transform it back by multiplying by X.
                Y = binop.mod(binop.multi(X,Y),self.mod)
                #Safety Check.
                LHS = binop.add(binop.mod(binop.square(Y),self.mod),binop.mod(binop.multi(X,Y),self.mod))
                if (LHS == RHS):
                    print("Point verified.")
                return [X,Y]
            else:
                nIterate = nIterate+1
                
    # Point operations
    #-------------------------------------------------------------
    def pointNeg(self,P):
        """Returns -P."""
        YRes = binop.add(P[0],P[1])
        return [P[0],YRes]

    def pointDbl(self,P):
        """Returns point 2P"""
        a = self.a
        theMod = self.mod
        if P == self.pointNeg(P):
            return O
        else:
            Lambda=binop.add(P[0],binop.divide(P[1],P[0],theMod))
            XRes=binop.add(binop.add(binop.square(Lambda),Lambda),a)
            XRes=binop.mod(XRes,theMod)
            YRes=binop.add(binop.add(binop.square(P[0]),binop.multi(Lambda,XRes)),XRes)
            YRes=binop.mod(YRes,theMod)
            return [XRes,YRes]

    def pointAdd(self,P,Q):
        """Returns point P+Q"""
        a = self.a
        theMod = self.mod
        if P != Q:
            if P == self.pointNeg(Q):
                return O
            elif P == O:
                return Q
            elif Q == O:
                return P
            else:
                Lambda=binop.divide(binop.add(P[1],Q[1]),binop.add(P[0],Q[0]),theMod)
                XRes=binop.add(binop.add(binop.add(binop.add(binop.square(Lambda),Lambda),P[0]),Q[0]),a)
                XRes=binop.mod(XRes,theMod)
                YRes=binop.add(binop.add(binop.multi(binop.add(P[0],XRes),Lambda),XRes),P[1])
                YRes=binop.mod(YRes,theMod)
                return [XRes,YRes]
        else:
            return self.pointDbl(P)

    def pointMul(self,P,n):
        """Returns point nP"""
        #Uses the double-and-add method for efficient point multiplication
        #Given a point P, what is nP?
        #Find the binary form of n
        bin_n = bin(n)[2:]
        adder = P
        result = O
        for i in range(len(bin_n)):
            if bin_n[len(bin_n)-i-1]=='1':
                result = self.pointAdd(adder,result)
            else:
                result = result
            adder = self.pointDbl(adder)
        return result
        
    def isNeg(P,Q):
        """Checks whether P = -Q"""
        if Q == self.pointNeg(P):
            return True
        else:
            return False
        
    # Calculating point order
    #-------------------------------------------------------------
    def pointorder(self,P):
        print("Calculating order of",P,"...")
        #Initialize the Hasse interval
        hasseL = math.ceil(fastexp(math.sqrt(fastexp(2,self.deg))-1,2))
        hasseU = math.ceil(fastexp(math.sqrt(fastexp(2,self.deg))+1,2))
        #From M in the Hasse interval, search M such that MP = O.
        for i in range(hasseL,hasseU):
            Q = self.pointMul(P,i)
            if Q == O:
                M = i
                break
            else:
                pass
        #Search whether a factor of M also annihilates P.
        for j in Primes:
            if j > M:
                break
            else:
                while M%j==0:
                    temp = M//j
                    #  Check whether M/j annihilates P.
                    ## If so, replace M with M/j. If this is not the case,
                    ## replace with another prime.
                    Q = self.pointMul(P,temp)
                    if Q == O:
                       M = temp
                    else:
                        break
        return M

#>>>>>>> Begin Simplified ECIES
# Setup: Parameters. (Public keys)
# The curve.
C = Curve(25)
print("=================================================")
print("a : ",C.a)
print("b : ",C.b)
print("=================================================")
P = C.generatePoint()
n = C.pointorder(P)
print("And the point is: ",P)
print("Its order is : ",n)
# Select two private keys. Calculate point Q.
m = random.randrange(1,n)
Q = C.pointMul(P,m)
k = random.randrange(1,n)
# Calculate kP and kQ.
kP = C.pointMul(P,k)
kQ = C.pointMul(Q,k)
# All set! Ask for a message.
plaintext = input("Input plaintext : ")
# Convert the plaintext into an array of ASCII codes in their binary form
plaintext_ascii = []
for text in plaintext:
    to_append = bin(ord(text))[2:]
    #Is the to_append's length is 7? If not, append zeroes in front.
    if (len(to_append) != 7):
        missing_digits = 7 - len(to_append)
        zeroes = ''
        for i in range(missing_digits):
            zeroes = zeroes + '0'
        to_append = zeroes + to_append
    plaintext_ascii.append(to_append)
# Separate the plaintext into blocks.
# The block length is calculated as follows: floor(GF degree / 7)
blocklength = C.deg // 7
plaintext_upper_limit = math.ceil(len(plaintext_ascii)/blocklength)
remainder = len(plaintext_ascii)%blocklength
plaintext_to_encrypt = []
for i in range(plaintext_upper_limit):
    to_append = ''
    if ((i != (plaintext_upper_limit - 1)) or (remainder == 0)):
        for j in range(blocklength):
            to_append = to_append + plaintext_ascii[blocklength*i + j]
    else:
        for j in range(remainder):
            to_append = to_append + plaintext_ascii[blocklength*i + j]            
    plaintext_to_encrypt.append(to_append)
# Mask every letter of the plaintext using the ECIES.
ciphertext = kP
for text in plaintext_to_encrypt:
    to_append = binop.mod(binop.multi(text,kQ[0]),C.mod)
    ciphertext.append(to_append)
print("Ciphertext : ",ciphertext)
# That's it for encryption! Now let's try decryption.
# Decrypt it using the private key.
plaintext_decrypt = []
sep = ""
my = C.pointMul(kP,m)
#print("my",my)
for i in range(2,len(ciphertext)):
    to_decode = binop.divide(ciphertext[i],my[0],C.mod)
    #If the string to be decoded is NOT divisible by 7, append zeroes to avoid errors.
    if (len(to_decode)%7 != 0):
        missing = 7 - (len(to_decode)%7)
        zeroes = ''
        for i in range(missing):
            zeroes = zeroes + '0'
        to_decode = zeroes + to_decode
    #Slice and dice.
    for i in range(len(to_decode)//7):
        start_slice = 7*i
        stop_slice = 7*(i+1)
        to_convert = int(to_decode[start_slice:stop_slice],2)
        to_append = chr(to_convert)
        plaintext_decrypt.append(to_append)
original = sep.join(plaintext_decrypt)
print("Decrypted plaintext : ",original)

