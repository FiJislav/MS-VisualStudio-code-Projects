cislo=input("Zadej cislo: ").strip()

if cislo.isdigit() == True:
    cislo = int(cislo)
    if cislo > 0:
        for i in range(0,100):
            print(i)
    elif cislo == 0:
        print("Cislo je 0")
    elif cislo < 0:
        for i in range(100, -1, -1):
            print(i)
    else:
        for i in range(100,0, -1):
            print(i)
else:
    print("neni cislo")


