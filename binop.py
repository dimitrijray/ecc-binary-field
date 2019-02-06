"""This module contains the necessary binary-field operations for applications in elliptic curve cryptography"""

def div(num,numR):
    """Returns the quotient of the polynomial division num/numR"""
    P=int(num,2)
    R=int(numR,2)
    degP=len(num)-1
    degR=len(numR)-1
    Q=0
    prevDegree=degP-degR
    for i in range(prevDegree,-1,-1):
        setDeg=degP-degR
        Q=Q<<1
        if i==setDeg:
            R_1=R<<setDeg
            P=P^R_1
            degP=len(bin(P)[2:])-1
            Q=Q^1
    return bin(Q)[2:]

def mod(num,numR):
    """Returns the remainder of the polynomial division num/numR"""
    P=int(num,2)
    R=int(numR,2)
    degP=len(num)-1
    degR=len(numR)-1
    if degR!=0:
        while (degR<=degP):
            setDeg=degP-degR
            R_1=R<<setDeg
            P=P^R_1
            degP=len(bin(P)[2:])-1
    else:
        P=0
    return bin(P)[2:]

def add(numA,numB):
    """Returns numA + numB"""
    #This algorithm will do bitwise-XOR-ing between A and B
    A=int(numA,2)
    B=int(numB,2)
    result = A^B
    return bin(result)[2:]

def multi(numA,numB):
    """Returns numA * numB"""
    A = int(numA,2)
    B = int(numB,2)
    #This algorithm will always multiply A with B.
    #Initialization step.
    if numB[len(numB)-1]=='0':
        C=0
    elif numB[len(numB)-1]=='1':
        C=A
    #Begin multiplication step.
    D=A
    for i in range(1,len(numB)):
        D=D<<1
        if numB[len(numB)-(i+1)]=='1':
            C=C^D
        else:
            C=C
    return bin(C)[2:]

def square(numA):
    """Returns numA^2"""
    #Quick squaring: just add zeroes.
    if numA[0]=='0':
        B=0
    elif numA[0]=='1':
        B=1
    for i in range(1,len(numA)):
        B=B<<2
        if numA[i]=='1':
            B=B^1
    return bin(B)[2:]

def euclid(numA,numB):
    """Returns gcd(numA,numB)"""
    #Initialization
    A=numA
    B=numB
    R=B
    while R!='0':
        R=mod(A,B)
        A=B
        B=R
    return A

def eeuclid(numA,numB):
    """Returns gcd(numA,numB), the coefficient of A, and the coefficient of B"""
    #Initialization
    A=numA
    B=numB
    R=B
    p1='1'
    p2='0'
    q1='0'
    q2='1'
    #Begin euclidean algorithm
    while R!='0':
        R=mod(A,B)
        Q=div(A,B)
        A=B
        B=R
        pA=add(p1,multi(Q,p2))
        p1=p2
        p2=pA
        qB=add(q1,multi(Q,q2))
        q1=q2
        q2=qB
    return [A,p1,q1]

def inverse(thePol,theMod):
    """Returns inverse of thePol modulo theMod"""
    #Initialization
    A=thePol
    B=theMod
    R=B
    p1='1'
    p2='0'
    #God rest ye, merry gentlemen.
    while R!='1':
        R=mod(A,B)
        Q=div(A,B)
        A=B
        B=R
        pA=add(p1,multi(Q,p2))
        p1=p2
        p2=pA
    return pA

def divide(numA,numB,theMod):
    """Returns numA * (numB)^-1 (mod theMod)"""
    #This algorithm will always divide A with B.
    #That is, binDivide(A,B,M) will return:
    #A * B(-1) (mod M)
    result = mod(multi(numA,inverse(numB,theMod)),theMod)
    return result

def binpow(num,exp):
    """Returns num^exp"""
    expo=bin(exp)[2:]
    mult=num
    result='1'
    for i in range(len(expo)-1,-1,-1):
        if expo[i]=='1':
            result=multi(mult,result)
        else:
            result=result
        mult=square(mult)
    return result

def powMod(num,exp,m):
    """Returns num^exp (mod m)"""
    expo=bin(exp)[2:]
    mult=num
    result='1'
    for i in range(len(expo)-1,-1,-1):
        if expo[i]=='1':
            result=multi(mult,result)
            result=mod(result,m)
        else:
            result=result
        mult=square(mult)
    return result

def isIrred(poly):
    """Checks whether polynomial poly is irreducible"""
    deg=len(poly)-1
    u='10'
    status=1
    for i in range(deg//2):
        u = mod(square(u),poly)
        d = euclid(poly,add(u,'10'))
        if d!='1':
            status=0
            break
    return status

def genIrred(n):
    """Generates an irreducible polynomial of degree n"""
    #Brute force trinomials.
    for i in range(1,n):
        thePoly=1
        NextDeg=n-i
        thePoly=(thePoly<<NextDeg)^1
        thePoly=(thePoly<<i)^1
        if isIrred(bin(thePoly)[2:]):
            return bin(thePoly)[2:]
            break
        else:
            pass
    for i in range(3,n):
        for j in range(2,i):
            for k in range(1,j):
                #polinom x^n+x^i+x^j+x^k+1
                thePoly=1
                NextDeg=n-i
                thePoly=(thePoly<<NextDeg)^1
                NextDeg=i-j
                thePoly=(thePoly<<NextDeg)^1
                NextDeg=j-k
                thePoly=(thePoly<<NextDeg)^1
                thePoly=(thePoly<<k)^1
                if isIrred(bin(thePoly)[2:]):
                    return bin(thePoly)[2:]
                    break
                else:
                    pass       
    return "N/A"
