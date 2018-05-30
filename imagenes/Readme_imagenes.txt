//GRAFICAS EXIGIDAS

Para prk4 se uso dt=T/6000, y para pteuler se uso dt=T/24000
en ambos caso se uso:
N=(tend-tini)/dt	//donde N es el numero de pasos
pos = {0.994, 0.0, 0.0}
vel = {0.0, -2.0015851063790825224053786224, 0.0}


Para pdp54 y dormandprince-points se uso:
dt=T/6000
N=134
hmin=0.0001
hmax=0.137
pos = {0.994, 0.0, 0.0}
vel = {0.0, -2.0015851063790825224053786224, 0.0}

Para plotgraf se usaron las condiciones anteriores


//GRAFICAS PARA VERIFICAR LA PERIODICIDAD DE LA ORBITA

Para periodicoeuler se uso: "se obtuvo 1 periodo completo"
dt=T/240000.0
N=(tend-tini)/dt
pos = {0.994, 0.0, 0.0}
vel = {0.0, -2.0015851063790825224053786224, 0.0}

Para periodicork4 se uso: "se obtuvieron 3 periodos completos"
dt=T/100000.0
N=(tend-tini)*3/dt
pos = {0.994, 0.0, 0.0}
vel = {0.0, -2.0015851063790825224053786224, 0.0}

Para periodicodp54 se uso: "se obtuvieron 6 periodos completos"
dt=0.0002
hmin=0.0000001
hmax=0.0002
N=400000
pos = {0.994, 0.0, 0.0}
vel = {0.0, -2.0015851063790825224053786224, 0.0}



//GRAFICAS CON Z DIFERENTE DE CERO

Para z01.1 y z01.2 se uso:
dt=0.0002
hmin=0.0000001
hmax=0.0002
N=150000
pos = {0.994, 0.0, 0.0}
vel = {0.5, -2.0015851063790825224053786224, 0.8}

Para z02.1 se uso:
dt=0.0002
hmin=0.0000001
hmax=0.0002
N=110000
pos = {0.994, 0.0, 0.0}
vel = {0.1, -2.0015851063790825224053786224, 0.1}

