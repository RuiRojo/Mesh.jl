using LinearAlgebra

export loadmesh_Lavia

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	Código auxiliar para geometría ( meshes y etc.)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# v. 13-8-2018 	
# Este código es común para cálculos de BEM y de KA

# Básicamente cargar meshes desde archivos de texto (STL y DAT) en diversos formatos
# Incorporé GMF que me permite salvar meshes con 64 bits (finalmente)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas para leer meshes de triángulos planos (STL texto ASCII) y salvar en formato VT 
#   	(vector de triángulos). Cada triángulo es una matriz de 3x4
#   	(optimización para meshes monstruosas)

# Función que levanta una mesh triangular plana 'file' STL y genera una estructura que es un vector
# de vectores 'MatrizTriangulos' de 'N x 4' donde cada fila contiene los cuatro vectores que
# caracterizan al triángulo (3 vértices y una normal). No existe la estructura usual de 'selv'
# y 'vertex' que es más económica en cuanto a espacio ocupado en memoria pero no en cuanto a
# rapidez de cálculo.
# 
function GetMesh_VT( file::String )
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Primera apertura para contar las líneas
    f = open(  file ) ;
    NL = countlines( f ) ;
    close( f ) ;
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Definición de parámetros y entes
    NE = Int64( (NL-2)/7 ); # Se descuentan las 'solid' y 'endsolid'
    verticeA = zeros(Float64,3) ;
    verticeB = zeros(Float64,3) ;
    verticeC = zeros(Float64,3) ;
    valorN = zeros(Float64,3) ;
    VecTriang = Array{ Array{Float64,2} }( NE ) ;
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Segunda apertura del archivo
    f = open( file ) ;
    if ( split( readline( f ) )[ 1 ] != "solid" ) # Leo 'solid'
        println( "Error: no existe línea 'solid'." ) ;
    end
    cara = 1 ;
    # Se identifican los vértices de cada elemento
    while !eof( f ) && cara <= NE # Se lee mientras no se llega al fondo
        subline = 1 ;
        while subline < 8
            linea = readline( f ) ; # Se lee una línea #line += 1 ;
            if subline == 1 # Normal del elemento
                valorN = parse.( Float64, split( linea )[ 3: end] ) ;
            end
            if subline == 3 # Vértice 1
                verticeA = parse.( Float64, split( linea )[ 2: end] ) ;
            end
            if subline == 4 # Vértice 2
                verticeB = parse.( Float64, split( linea )[ 2: end] ) ;
            end
            if subline == 5 # Vértice 3
                verticeC = parse.( Float64, split( linea )[ 2: end] ) ;
            end
            subline += 1 ;
        end # Armo la matriz con los vértices del triángulo y la normal
        VecTriang[ cara ] = hcat( verticeA, verticeB, verticeC, valorN ) ; 
        cara += 1 ;
    end
    close( f ) ; # Cierre del archivo
    return VecTriang
end

# Funcion para guardar un mesh tipo VT (vector de matrices) de manera binaria
function WriteMeshBin_VT( file::String, VT::AbstractArray )
    open(file, "w") do f_out
        for i in VT 
            for j in eachindex( i )
                write( f_out, i[ j ] ) ;
            end
        end
    end
end

# Función para leer una mesh de tipo VT (vector de matrices)
function ReadMeshBin_VT( file::String, dtype::Type )
    VT = Array{ Array{Float64, 2 }}( 0 ) ;
    open(file, "r") do f_in 
        while !eof( f_in ) 
            push!( VT, reshape( read( f_in, dtype, 12 ), 3, 4 )  ) ;
        end
    end
    return VT
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas para leer meshes binarias y generar un matriz de puntos 'MP'
#

# Función para leer una mesh de tipo MP (matriz de puntos)
# Resultó esto más óptimo que leer y asignar a un vector y luego aplicar
# un reshape al vector. No pude generar matriz directamente y asignar a ese ente
# Si bien se exporta una matriz, se debe prealocar un vector de vectores para ir
# leyendo secuencialmente de a un triángulo.
function ReadMeshBin_MP( file::String, dtype::Type )
    MP = Array{ Array{ Float64, 2 }}(undef,  0 ) ; 
    open(file, "r") do f_in 
        while !eof( f_in ) 
            push!( MP, read!(f_in, Array{dtype}(undef, 12, 1) )) ;
        end
    end
    return hcat( MP... ) ;
end

# Función que convierte las estructuras 'selv', 'vertex', 'normales' de meshes de triángulos planas
# en una matriz 'MP'
function SelVertNorm_to_MP( selv, vertex, normales )
    N = size(selv)[1] ;
    MP = Array{Float64}(12,N) ;
    for i = 1 : N
        MP[1:3,i] = vertex[selv[i,1],:] ;
        MP[4:6,i] = vertex[selv[i,2],:] ;
        MP[7:9,i] = vertex[selv[i,3],:] ;
        MP[10:12,i] = normales[i,:] ;
    end
    return MP
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas para trabajar con meshes binarias salvadas con FileIO y JLD
#

# Se encontró que no es óptimo trabajar con estructuras que no sean de 'bits type' (matrices, vectores)
# porque al guardar-leer binariamente estructuras complejas (del tipo vector de vectores) se hace con
# un overhead intolerable puesto que se almacenan los punteros (o algo así --ask to Rui--).

# Función para guardar binariamente en archivo 'file' las estructuras de un mesh.
# Requieren packages 'FileIO' y 'JLD'
function savemesh( file::String, norms::Array, selv::Array, vert::Array )
    save( file, "normales", norms, "selv", selv, "vertex", vert ) ;
end

# Función para recuperar las estructuras de un 'mesh' salvadas en archivo.
# Requieren packages 'FileIO' y 'JLD'
function LevantarMesh( file::String )
    return normales, selv, vert  = load( file, "normales", "selv", "vertex" ) ;
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas para leer meshes cuadráticas (DAT texto ASCII) y salvar en estructuras 'selv' 'vertex'
#

# Función que levanta una mesh cuadrática (*.dat) y devuelve lista de elementos y vértices
# como matrices. Que la salida sea una matriz en lugar de un vector de vectores es importante 
# a la hora de convertirlo a shared array.
# input:
#	file		:	Archivo con el mesh
# output:
#	selv		:	Lista de elementos (matriz de N x 6 )
#	vert		:	Lista de vértices (matriz de N x 3 )
function GetMeshQuad( file::String )
    f = open(  file ) ;
    F = readlines( f ) ;
    SizeF = size(F)[1] ;
    close( f ) ;
    NroVertices = parse.( Int64, split( F[1] ) )[1] ;
    NroElementos = parse.( Int64, split( F[1] ) )[2] ;
    VerticesLista = Array{Float64}(undef, NroVertices, 3 ) ;
    Selv = zeros( Int64, NroElementos, 6 ) ;
    contador = Int64( 0 ) ;
    for i = 2 : NroVertices + 1
        VerticesLista[i-1,:] = transpose( parse.( Float64, split( F[i] )[2:end] ) ) ;
    end
    for i = NroElementos + 2 : SizeF
        Linea = parse.( Float64, split( F[i] ) ) ;
        # if Linea[2] == 103
        # 	;
        if Linea[2] == 206 # La primera vez que caíste acá
            contador += 1;
            Selv[ contador,: ] = convert.( Int64, Linea[3:end] ) ;
        end 
    end
    println( "Nro. de triángulos : ", contador )
    NroFaces = contador ;
    Selv2 = Selv[1:NroFaces,:] ;
    close( f )
    return Selv2, VerticesLista
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas para graficar meshes en GNUPLOT
#

# Función que salva una mesh para su graficación en gnuplot salvando separadamente
# los segmentos (c/u de los cuales son líneas consecutivas). Cada triángulo de vértices
# A,B,C se guarda con el esquema :
# 	A_x A_y A_z
#	B_x B_y B_z
#	
# 	A_x A_y A_z
#	C_x C_y C_z
#
# 	B_x B_y B_z
#	C_x C_y C_z
#	
function MeshForPlot( file::String, Selv::AbstractArray, Vert::AbstractArray )
    open(file, "w") do f_out # Abro el archivo
        for i = 1 : size( Selv )[1] # Recorro selv
            writedlm( f_out, [ Vert[Selv[i,1],:] ] ) ;
            writedlm( f_out, [ Vert[Selv[i,2],:] ] ) ;
            writedlm( f_out, '\n') ;
            writedlm( f_out, [ Vert[Selv[i,1],:] ] ) ;
            writedlm( f_out, [ Vert[Selv[i,3],:] ] ) ;
            writedlm( f_out, '\n') ;
            writedlm( f_out, [ Vert[Selv[i,2],:] ] ) ;
            writedlm( f_out, [ Vert[Selv[i,3],:] ] ) ;
            writedlm( f_out, '\n') ;
        end
    end
end

# Se salvarán los vértices y un punto intermedio extra
function VerticesForPlotInterp( file::String, Selv::AbstractArray, Vert::AbstractArray )
    open(file, "w") do f_out # Abro el archivo
        for i = 1 : size( Selv )[1] # Recorro triángulos (selv)
            writedlm( f_out, [ Vert[Selv[i,1],:] ] ) ;
            writedlm( f_out, [ Vert[Selv[i,2],:] ] ) ;
            writedlm( f_out, [ Vert[Selv[i,3],:] ] ) ;
            writedlm( f_out, [ 1/3 * ( Vert[Selv[i,1],:] + Vert[Selv[i,2],:] + Vert[Selv[i,3],:] ) ] ) ;
            writedlm( f_out, [ Vert[Selv[i,4],:] ] ) ;
            writedlm( f_out, [ Vert[Selv[i,5],:] ] ) ;
            writedlm( f_out, [ Vert[Selv[i,6],:] ] ) ;
        end
    end
end

# Se salvarán los vértices y un punto intermedio extra
# falta encapsular las funciones para que sea más operativo esto.
function VerticesForPlotInterp_Best( file::String, Selv::AbstractArray, Vert::AbstractArray )
    Values = t -> real.( exp.( im * 5 * t[1] ) ) ;
    Prom3 = ( a,b,c ) -> 1/3 * ( a + b + c ) ;
    Prom2 = ( a,b ) -> 1/2 * ( a + b ) ;
    Prom2V1 = ( a,b ) -> ( 1/6 * a + 5/6 * b ) ;
    Prom2V2 = ( a,b ) -> ( 2/3 * a + 1/3 * b ) ;
    open(file, "w") do f_out # Abro el archivo
        for i = 1 : size( Selv )[1] # Recorro triángulos (selv)
            V1 = Vert[Selv[i,1],:] ;
            V2 = Vert[Selv[i,2],:] ;
            V3 = Vert[Selv[i,3],:] ;
            V4 = Vert[Selv[i,4],:] ;
            V5 = Vert[Selv[i,5],:] ;
            V6 = Vert[Selv[i,6],:] ;
            writedlm( f_out, hcat( V1... , Values( V1 ) ) ) ;
            writedlm( f_out, hcat( V2... , Values( V2 ) ) ) ;
            writedlm( f_out, hcat( V3... , Values( V3 ) ) ) ;
            writedlm( f_out, hcat( Prom3( V1, V2, V3)... , Values( Prom3( V1, V2, V3) ) ) ) ;
            writedlm( f_out, hcat( Prom2( V1, V5 )..., Values( Prom2( V1, V5 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2( V2, V6 )..., Values( Prom2( V2, V6 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2( V3, V4 )..., Values( Prom2( V3, V4 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2( V4, V5 )..., Values( Prom2( V4, V5 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2( V4, V6 )..., Values( Prom2( V4, V6 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2( V5, V6 )..., Values( Prom2( V5, V6 ) ) ) ) ;
            #	
            writedlm( f_out, hcat( Prom2V1( V1, V5 )..., Values( Prom2V1( V1, V5 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V1( V2, V6 )..., Values( Prom2V1( V2, V6 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V1( V3, V4 )..., Values( Prom2V1( V3, V4 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V1( V4, V5 )..., Values( Prom2V1( V4, V5 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V1( V4, V6 )..., Values( Prom2V1( V4, V6 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V1( V5, V6 )..., Values( Prom2V1( V5, V6 ) ) ) ) ;
            #	
            writedlm( f_out, hcat( Prom2V2( V1, V5 )..., Values( Prom2V2( V1, V5 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V2( V2, V6 )..., Values( Prom2V2( V2, V6 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V2( V3, V4 )..., Values( Prom2V2( V3, V4 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V2( V4, V5 )..., Values( Prom2V2( V4, V5 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V2( V4, V6 )..., Values( Prom2V2( V4, V6 ) ) ) ) ;
            writedlm( f_out, hcat( Prom2V2( V5, V6 )..., Values( Prom2V2( V5, V6 ) ) ) ) ;
            writedlm( f_out, hcat( V4... , Values( V4 ) ) ) ;
            writedlm( f_out, hcat( V5... , Values( V5 ) ) ) ;
            writedlm( f_out, hcat( V6... , Values( V6 ) ) ) ;
        end
    end
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas para trabajar con los seis parches de la esfera
#

# Función que devuelve 'x,y,z' viviendo en la esfera dado
# el valor del parche 'pat' y las coordenadas 'u', 'v' de
# parametrización. 
#		-1 <= u, v <= 1
function Parche( pat::Int, u, v )
    fact = sqrt.( 1 .+ u.^2 + v.^2 ) ;
    if pat == 1
        return	u./fact, v./fact, 1 ./fact ;
    elseif pat == 2
        return	u./fact, 1 ./fact, v./fact ;
    elseif pat == 3
        return	1 ./fact, v./fact, u./fact ;
    elseif pat == 4
        return	u./fact, v./fact, -1 ./fact ;
    elseif pat == 5
        return	u./fact, -1 ./fact, v./fact ;
    elseif pat == 6
        return	-1 ./fact, v./fact, u./fact ;
    end
end

# Función que devuelve una matriz de vectores de incidencia a partir de un
# set de vectores de discretización -1 <= u, v <= 1 para el parche 'patch'
# entre 1 y 6.
function MatIncParche( patch, u, v )
    M = size( u )[ 1 ] ;
    N = size( v )[ 1 ] ;
    Angulos = Array{ Tuple{ Float64, Float64 }, 2 }( M, N ) ;
    for i = 1 : M
        for j = 1 : N
            (x,y,z) = Parche( patch, u[i], v[j] ) ;
            phi = atan2( y, x ) ;
            theta = atan2( sqrt( x^2 + y^2 ), z ) ;
            Angulos[ i, j ] = (theta, phi) ;
        end
    end
    return GetKinc.( Angulos ) ;
end

function MatIncParche2( patch, u, v )
    M = size( u )[ 1 ] ;
    N = size( v )[ 1 ] ;
    Angulos = Array{ Tuple{ Float64, Float64 }, 2 }( M, N ) ;
    for i = 1 : M
        for j = 1 : N
            (x,y,z) = Parche( patch, u[i], v[j] ) ;
            phi = atan2( y, x ) ;
            theta = atan2( sqrt( x^2 + y^2 ), z ) ;
            Angulos[ i, j ] = (theta, phi) ;
        end
    end
    return ( Angulos ) ;
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas para leer meshes de triángulos planos (STL texto ASCII) y salvar en estructuras 'selv' 'vertex'
#

# Función que levanta una mesh lineal STL en 'file' y devuelve las estructuras
#	normales, selv, vertex
# Versión acelerada para operar en meshes grandes pero no monstruosos.
# En meshes monstruosos se opta por utilizar otro tipo de estructuras.
function GetMesh_fast( file::String )
    # Lectura del archivo .stl #
    if file[end-3:end] != ".stl"
        println("Error : No parecece ser un archivo STL") ;
    end
    line = 2
    valorN = zeros(Float64,1,3) ;
    # Primera apertura para contar las líneas
    f = open(  file ) ;
    NL = countlines( f ) ;
    close( f ) ;
    NE = Int64( (NL-2)/7 ); # Se descuentan las 'solid' y 'endsolid'
    indicesCaras = zeros( Int64, NE, 3 ) ;
    Normales = zeros( Float64, NE, 3 ) ;
    valorA = zeros(Float64,1,3) ;
    valorB = zeros(Float64,1,3) ;
    valorC = zeros(Float64,1,3) ;
    nro_vert_A = Int64( 0 ) ;
    nro_vert_B = Int64( 0 ) ;
    nro_vert_C = Int64( 0 ) ;
    VerticesLista = Array{Array{Float64,2}}(undef, 0) ;

    # Segunda apertura del archivo para construir la lista de vértices y normales
    f = open( file ) ;
    if ( split( readline( f ) )[ 1 ] != "solid" ) # Leo 'solid'
        println( "Error: no existe línea 'solid'." ) ;
    end
    cara = 1 ;
    # Se llena la lista de vértices sin chequear los repetidos
    while !eof( f ) && cara <= NE # Se lee mientras no se llega al fondo
        subline = 1 ;
        while subline < 8
            linea = readline( f ) ; # Se lee una línea #line += 1 ;
            if subline == 1 # Normal del elemento
                valorN = hcat( parse.( Float64, split( linea )[ 3: end] )... ) ;
            end
            if subline == 3 # Vértice 1
                valorA = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                push!( VerticesLista, valorA ) ;
            end
            if subline == 4 # Vértice 2
                valorB = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                push!( VerticesLista, valorB ) ;
            end
            if subline == 5 # Vértice 3
                valorC = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                push!( VerticesLista, valorC ) ;
            end
            subline += 1 ;
        end
        Normales[ cara, : ] = valorN; # Normales de los elementos
        cara += 1 ;
    end
    close( f ) ;
    # Cierre del archivo
    println("Primer cierre") ;

    VerticesListaClean = CleanRepetidos( VerticesLista ) ; # Limpieza de repetidos
    VLCcolumn1 = vcat(VerticesListaClean...)[:,1] ; # Se transforma a matriz y se toma la primer columna
    println("Limpieza repetidos") ;
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Tercera apertura del archivo
    f = open( file ) ;
    if ( split( readline( f ) )[ 1 ] != "solid" ) # Leo 'solid'
        println( "Error: no existe línea 'solid'." ) ;
    end
    cara = 1 ;
    # Se identifican los vértices de cada elemento
    while !eof( f ) && cara <= NE # Se lee mientras no se llega al fondo
        subline = 1 ;
        while subline < 8
            linea = readline( f ) ; # Se lee una línea #line += 1 ;
            if subline == 3 # Vértice 1
                valorA = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                nro_vert_A = GetIdxVertex( VerticesListaClean, VLCcolumn1, valorA ) ;
            end
            if subline == 4 # Vértice 2
                valorB = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                nro_vert_B = GetIdxVertex( VerticesListaClean, VLCcolumn1, valorB ) ;
            end
            if subline == 5 # Vértice 3
                valorC = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                nro_vert_C = GetIdxVertex( VerticesListaClean, VLCcolumn1, valorC ) ;
            end
            subline += 1 ;
        end
        indicesCaras[ cara, : ] = [ nro_vert_A, nro_vert_B, nro_vert_C ];
        cara += 1 ;
        if mod( cara, 10000 ) == 0
            println("Caras :", cara )
        end
    end
    close( f ) ;
    # Cierre del archivo
    println("Número de triángulos / vértices: ", size(indicesCaras)[1]," / ",size(VerticesListaClean)[1] ) ;
    # Falta la parte de chequeo de orientación y normales
    return Normales, indicesCaras, vcat( VerticesListaClean... )
end
function GetMesh_faster( file::String )
    println("FASTER yeah")
    # Lectura del archivo .stl #
    if file[end-3:end] != ".stl"
        println("Error : No parecece ser un archivo STL") ;
    end
    line = 2
    valorN = zeros(Float64,1,3) ;
    # Primera apertura para contar las líneas
    lines = readlines(file)
    NL = length(lines)
    NE = Int64( (NL-2)/7 ); # Se descuentan las 'solid' y 'endsolid'
    indicesCaras = zeros( Int64, NE, 3 ) ;
    Normales = zeros( Float64, NE, 3 ) ;
    valorA = zeros(Float64,1,3) ;
    valorB = zeros(Float64,1,3) ;
    valorC = zeros(Float64,1,3) ;
    nro_vert_A = Int64( 0 ) ;
    nro_vert_B = Int64( 0 ) ;
    nro_vert_C = Int64( 0 ) ;
    VerticesLista = Array{Array{Float64,2}}(undef, 0) ;

    # Segunda apertura del archivo para construir la lista de vértices y normales
    cnt = 0
    if ( split( lines[cnt+=1] )[ 1 ] != "solid" ) # Leo 'solid'
        println( "Error: no existe línea 'solid'." ) ;
    end
    cara = 1 ;
    # Se llena la lista de vértices sin chequear los repetidos
    while cnt<NL && cara <= NE # Se lee mientras no se llega al fondo
        subline = 1 ;
        while subline < 8
            linea = lines[cnt+=1]; # Se lee una línea #line += 1 ;
            if subline == 1 # Normal del elemento
                valorN = hcat( parse.( Float64, split( linea )[ 3: end] )... ) ;
            end
            if subline == 3 # Vértice 1
                valorA = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                push!( VerticesLista, valorA ) ;
            end
            if subline == 4 # Vértice 2
                valorB = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                push!( VerticesLista, valorB ) ;
            end
            if subline == 5 # Vértice 3
                valorC = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                push!( VerticesLista, valorC ) ;
            end
            subline += 1 ;
        end
        Normales[ cara, : ] = valorN; # Normales de los elementos
        cara += 1 ;
    end
    # Cierre del archivo
    println("Primer cierre") ;

    VerticesListaClean = CleanRepetidos( VerticesLista ) ; # Limpieza de repetidos
    VLCcolumn1 = vcat(VerticesListaClean...)[:,1] ; # Se transforma a matriz y se toma la primer columna

    mat2tup(x) = (x[1], x[2], x[3])
    verticesVec = map(mat2tup, VerticesListaClean)
    vert2Idx = Dict(zip(verticesVec, 1:length(verticesVec)))
    getVertIdx(vert::Matrix) = getVertIdx(mat2tup(vert))
    getVertIdx(vert::Tuple) = vert2Idx[vert]

    println("Limpieza repetidos") ;
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Tercera apertura del archivo
    cnt=0
    if ( split( lines[cnt+=1] )[ 1 ] != "solid" ) # Leo 'solid'
        println( "Error: no existe línea 'solid'." ) ;
    end
    cara = 1 ;
    # Se identifican los vértices de cada elemento
    while cnt < NL && cara <= NE # Se lee mientras no se llega al fondo
        subline = 1 ;
        while subline < 8
            linea = lines[cnt+=1] ; # Se lee una línea #line += 1 ;
            if subline == 3 # Vértice 1
                valorA = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                nro_vert_A = getVertIdx( valorA ) ;
            end
            if subline == 4 # Vértice 2
                valorB = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                nro_vert_B = getVertIdx( valorB ) ;
            end
            if subline == 5 # Vértice 3
                valorC = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                nro_vert_C = getVertIdx( valorC ) ;
            end
            subline += 1 ;
        end
        indicesCaras[ cara, : ] = [ nro_vert_A, nro_vert_B, nro_vert_C ];
        cara += 1 ;
        if mod( cara, 100000 ) == 0
            println("Caras :", cara )
        end
    end
 
    # Cierre del archivo
    println("Número de triángulos / vértices: ", size(indicesCaras)[1]," / ",size(VerticesListaClean)[1] ) ;
    # Falta la parte de chequeo de orientación y normales
    return Normales, indicesCaras, vcat( VerticesListaClean... )
end


# Función que toma una lista de vértices VertList (con vértices repetidos)
# y la 'limpia' exportando una lista de vértices no repetidos.
# 
# Como los vértices se pasan ordenados, los repetidos están todos juntos.
# Dado un vértice de esta lista ordenada 'VL', el único que puede ser igual
# es el anterior. Esto simplifica la búsqueda de repetidos notablemente.
function CleanRepetidos( VertList::AbstractArray )
    epsil = 1E-12 ; # Hardcoded precision (podríamos poner 0 puesto que provienen del mismo parseo)
    VL = sortslices( vcat(VertList...), dims=1 ) ; # Transformo a matriz y ordeno por la coordenada 'x'
    VLC = Array{ Array{ Float64,2 } }(undef, 0 ) ; # Lista vacía de vértices 'clean'
    push!( VLC, hcat( VL[1,:] )' ) ; # El primero entra siempre
    for i = 2 : size( VL )[1] # Recorro la lista
        if norm( VL[i,:] - VL[i-1,:] ) < epsil # Si fue igual al anterior es un repetido. No entra
            ; # Ya está en la lista
        else # Si no fue igual al anterior es la primera (tal vez única) ocurrencia. Entra
            push!( VLC, hcat( VL[i,:] )' ) ; 
        end 
    end 
    return VLC
end

# Primitiva manera de identificar un vértice 'V' en la lista 'VertList'. Chequea toda la lista y no 
# es por ello muy funcional.
function GetIdxVertexF( VertList::AbstractArray, V::AbstractArray )
    return find( q ->  norm( V - q ) == 0 , VertList )[ 1 ]
end

# Función que devuelve el 'id' numérico dentro de la lista de vértices 'VertList'
# de un dado vértice 'V' si se le pasa como entrada el vector 'VertList1'que 
# corresponde a la primer columna de 'VertList'.
# 'VertList' está ordenada por la primer componente de los vértices (la componente x).
function GetIdxVertex( VertList::AbstractArray, VertList1::AbstractArray, V::AbstractArray )
    # const X = vcat(VertList...)[:,1] ; # Vector que tiene las componentes 'x' (están ordenadas de menor a mayor)
    v1 = V[1] ; # Componente 'x' de V
    L = length( VertList1 ) ;
    idx = 1 ;
    while v1 > VertList1[idx] # Busco la primera ocurrencia
        idx += 1 ;
    end
    first = idx ;
    final = idx ; # Por defecto
    while ( final < L && v1 == VertList1[final+1] )
        final +=1 ;
    end
    for j = first:final 
        if norm( V - VertList[j]) == 0
            return j
        end 
    end 
    return 0 # Si ocurrió error retorno 0
end


# Función que levanta una mesh triangular plana. Versión tomada del código BEM
#
# Código original de Juan (EsferaOrientadaTipoKirkup.jl) con upgrades de Lavia
# para que tome archivos 'stl' más generales
# Ha sido supersedida por 'GetMesh_fast' puesto que el tiempo para cargar meshes
# crece con 'N' de manera geométrica dado el algoritmo de búsqueda de repetidos 
# utilizado.
function GetMesh( file::String ) 
    # Lectura del archivo .stl #
    line = 2
    normalStr = 0;
    valorN = zeros(Float64,1,3) ;
    # Apertura para contar las líneas
    f = open(  file ) ;
    NL = countlines( f ) ;
    close( f ) ;

    verticesInicial = zeros( NL, 3 ) ;
    NE = Int64( (NL-2)/7 ); # Se descuentan las 'solid' y 'endsolid'
    indicesCaras = zeros( Int64, NE, 3 ) ;
    Normales = zeros( Float64, NE, 3 ) ;
    valorA = zeros(Float64,1,3);
    valorB = zeros(Float64,1,3);
    valorC = zeros(Float64,1,3);
    nro_vert_A = 0 ;
    nro_vert_B = 0 ;
    nro_vert_C = 0 ;

    VerticesLista = Array{Array{Float64,2}}(0) ;

    # Nueva apertura del archivo
    f = open( file ) ;
    if ( split( readline( f ) )[ 1 ] != "solid" ) # Leo 'solid'
        println( "Error: no existe línea 'solid'." ) ;
    end 
    cara = 1 ;
    while !eof( f ) && cara <= NE # Se lee mientras no se llega al fondo
        subline = 1 ;
        while subline < 8
            linea = readline( f ) ; # Se lee una línea #line += 1 ;
            if subline == 1 # Normal del elemento
                valorN = hcat( parse.( Float64, split( linea )[ 3: end] )... ) ;
            end
            if subline == 3 # Vértice 1
                valorA = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                # Cargar el vértice en la lista e identificar su id
                nro_vert_A = VectorEnLista( VerticesLista, valorA ) ;
            end
            if subline == 4 # Vértice 2
                valorB = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                nro_vert_B = VectorEnLista( VerticesLista, valorB ) ; # Ver si está en la lista
            end
            if subline == 5 # Vértice 3
                valorC = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
                nro_vert_C = VectorEnLista( VerticesLista, valorC ) ;# Ver si está en la lista
            end
            subline += 1 ;
        end
        indicesCaras[ cara, : ] = [ nro_vert_A, nro_vert_B, nro_vert_C ];
        Normales[ cara, : ] = valorN;
#			centroidesCaras[cara,:]=(valorA+valorB+valorC)/3;
        cara += 1 ;
    end
    close( f ) ;
    # Cierre del archivo

    # Fin de lectura archivo STL

    ### Acomodo los Indices tipo Kirkup y chequeo las normales del stl
    indicesCarasOrientado = zeros( Int64, size(indicesCaras) )
    contador = 0;
    for i = 1 : NE
        A = VerticesLista[ indicesCaras[i,1], 1 ];
        B = VerticesLista[ indicesCaras[i,2], 1 ]; 
        C = VerticesLista[ indicesCaras[i,3], 1 ];
        n1 = cross(vec(B-A),vec(C-A));
        n1 = n1/norm(n1);
        centroide = ( A + B + C )/3 ;
        # sólo sirve para cuerpos convexos centrados (o algo así?)
            # if dot( vec( centroide ) , vec(Normales[i,:] ) ) < 0
            # correccion_normales[ i ] = 1 ; # Si tienen sentido diferente
            # Normales[ i, : ] = -Normales[i,:] ;
        #end
        # Comparo normal con el cálculo a partir de los vértices
        if ( dot( vec( n1 ), vec( Normales[i,:] ) ) > 0.99 ) # Paralelos
            indicesCarasOrientado[i,:] = indicesCaras[i,:] ;
            contador += 1 ;
        end
        if ( dot( vec( n1 ), vec( Normales[i,:] ) ) < -0.99 ) # Antiparalelos
            # Corrijo el orden de los vértices
            indicesCarasOrientado[i,1] = indicesCaras[i,1];
            indicesCarasOrientado[i,3] = indicesCaras[i,2];
            indicesCarasOrientado[i,2] = indicesCaras[i,3];
            # Corrijo el sentido de la normal
            Normales[ i, : ] = -Normales[i,:] ;
            println("Invertido : ", i);
            contador += 1 ;
        end
    end
    # si esta bien contador tiene que ser igual a numero de elementos
    println( ( contador == NE ) ) ;
    # Chequeo orientación correcta. La normal en 'Normales' coincide
    # con la normal calculada a partir de los vértices
    contador = 0
    for i = 1:NE
        A = VerticesLista[indicesCarasOrientado[i,1],1] ;
        B = VerticesLista[indicesCarasOrientado[i,2],1] ;
        C = VerticesLista[indicesCarasOrientado[i,3],1] ;
        n1 = cross(vec(B-A),vec(C-A)) ;
        n1 = n1/norm(n1) ;
        if ( abs( dot(vec(n1), vec(Normales[i,:]))-1 ) < 1e-3 )
            contador += 1 ;
        end
    end
    # Si está todo re bien, la siguiente variable es true.
    println("Bien ordenado: ", ( contador == NE ), " Ordenados/Elements :", contador, "/", NE ) ;
    return Normales, indicesCarasOrientado, vcat( VerticesLista... )
end

# Función para incorporar un vector de vértices 'V' a una lista 'VertList'
# si no está y devolver un id (su nro en la lista). Si V existe en VertList
# se devuelve su ubicación.
function VectorEnLista( VertList::AbstractArray, V::Array{Float64,2} )
    found = false ;
    idx = 0 ;
    eps = 1e-10 ;
    for i = 1 : length( VertList )
        if norm( V - VertList[i] ) < eps
            found = true ;
            idx = i ;
        end
    end
    if !found # Si no está se incorpora
        push!( VertList, V ) ;
        return length( VertList )
    else # Si ya estaba en la lista de vértices
        return idx
    end
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas para leer meshes de triángulos planos (texto ASCII) que utilizan el paquete 'Meshes'
#

# Se utiliza el paquete 'Meshes' para obtener las estructuras usuales necesarias
function GetMesh_Meshes( File::String )
    mesh = load( File ) ; # Carga del mesh

    SELV = vcat( map(x -> ( convert.( Int64, x ) )' , mesh.faces )... )
    VERTEX = vcat( map(x -> ( convert.( Float64, x ) )' , mesh.vertices )... )
    NORMALES = vcat( map(x -> ( convert.( Float64, x ) )' , mesh.normals )... )

    return NORMALES, SELV, VERTEX
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas para leer meshes de triángulos planos (texto ASCII) en formato GMF (sin normales)
#

# Las normales se calculan con otra rutina para tener una salida estandarizada con respecto
# a las rutinas para formato "STL"
function GetMesh_GMF( file::String )
    # Lectura del archivo .mesh #	
    if file[end-4:end] != ".mesh"
        println("Error : No parecece ser un archivo GMF (extensión '.mesh')") ;
    end
    VerticesLista = Array{Array{Float64,2}}(undef, 0) ;
    TrianglesLista = Array{Array{Int64,2}}(undef, 0) ;
    # Apertura del archivo para construir la lista de vértices y triángulos
    f = open( file ) ;
    # Se llena la lista de vértices sin chequear los repetidos
    while !eof( f ) # Se lee mientras no se llega al fondo
        linea = readline( f ) ; # Se lee una línea
        if linea == "Vertices" # Buscar los vértices
            nro_vertices = parse( Int64, readline( f ) ) ; # Número de vértices
            for i = 1 : nro_vertices
                vertex = hcat( parse.( Float64, split( readline( f ) )[ 1:3 ] )... ) ;
                push!( VerticesLista, vertex ) ;
            end
            if readline( f ) != ""
                println("ERROR : hay más vértices de los que debiera.")  ;
            end
        end
        if linea == "Triangles"
            nro_triangles = parse( Int64, readline( f ) ) ; # Número de triángulos
            for i = 1 : nro_triangles
                triangle = hcat( parse.( Int64, split( readline( f ) )[ 1:3 ] )... ) ;
                push!( TrianglesLista, triangle ) ;
            end
            if readline( f ) != ""
                println("ERROR : hay más triángulos de los que debiera.")  ;
            end
        end
    end
    close( f ) ;
    println("Número de triángulos / vértices: ", size( TrianglesLista )[1]," / ", size( VerticesLista)[1] ) ;
    Tri = vcat( TrianglesLista... ) ;
    Ver = vcat( VerticesLista... ) ;
    Nor = BuildNormales( Tri, Ver ) ;
    return Tri, Ver, Nor
end

# Función para construir las normales a partir de una lista de triángulos 'Tri' y los vértices 'Ver'
function BuildNormales( Tri, Ver )
    N = size(Tri)[1] ;
    Normales = Array{Array{Float64,2}}(undef, N ) ;
    for q = 1 : N
        Normales[ q ] = hcat( normal3( Ver[ Tri[q,1],:], Ver[Tri[q,2],:], Ver[Tri[q,3],:] )... ) ;
    end
    return vcat( Normales... )
end


# Función que da la normal a tres vectores. Es la normal externa
# porque los vectores se esperan en CCW order. Traducción rutina de S. Kirkup
function normal3( VECA, VECB, VECC )
    VBMVA = VECB - VECA ; # CALL SUBV3(VECB,VECA,VBMVA)
    VCMVA = VECC - VECA ; # CALL SUBV3(VECC,VECA,VCMVA)     
    CR1 = VBMVA[2]*VCMVA[3] - VBMVA[3]*VCMVA[2] ;
    CR2 = VCMVA[1]*VBMVA[3] - VBMVA[1]*VCMVA[3] ;
    CR3 = VBMVA[1]*VCMVA[2] - VBMVA[2]*VCMVA[1] ;
    ACR = sqrt( CR1*CR1 + CR2*CR2 + CR3*CR3 ) ;
    return [ CR1/ACR, CR2/ACR, CR3/ACR ] ;
end


##---------------------
# Esto ya es propio, para convertir los outputs de Ed en `<:AbstractMesh`



"""
	loadmesh_Lavia(fn, s)

`s` es un símbolo `:STL_bin`, `:STL_ASCII`, `:GMF` o `:Curvo`.

Equivalentemente, un string con la extensión à-la-Lavia: "bin", "stl", "mesh" o "dat".

Devuelve un `AbstractMesh`. O sea, type unstable.
"""
function loadmesh_Lavia(fn, s)
	d = Dict(
		"bin"  => :STL_bin,
		"stl"  => :STL_ASCII,
		"mesh" => :GMF,
		"dat"  => :Curvo	
	)

	s = lowercase(s)
	sym = d[s]

	if sym == :STL_bin
		loadmesh_Lavia_stl_bin(fn)

	elseif sym == :STL_ASCII
		loadmesh_Lavia_stl_ascii(fn)

	elseif sym == :GMF
		loadmesh_Lavia_gmf(fn)

	elseif sym == :Curvo
		loadmesh_Lavia_curvo(fn)
	end
end


function loadmesh_Lavia_stl_bin(fn) :: LightMesh
    s = ReadMeshBin_MP(fn, Float64)
    vs1, vs2, vs3 = s[[1; 4; 7], :], s[[2; 5; 8], :], s[[3; 6; 9],:]
    vss = map(Vertex, vs1, vs2, vs3)
    vss = mapslices(NTuple{3}, vss; dims=1)[1, :]

    return convert(LightMesh, SimpleMesh(vss))
end

faster() = true #prov, debería fijarse en true y ya

function loadmesh_Lavia_stl_ascii(fn) :: LightMesh
    norm, selv, vertex  = (faster() ? GetMesh_faster : GetMesh_fast)(fn)
    vs = Vertex.(mapslices(NTuple{3, Float32}, vertex; dims=2)[:, 1])
    faces = mapslices(NTuple{3, FaceIndex}, selv, dims=2)[:, 1]
     
    return LightMesh(vs, faces)
end

function loadmesh_Lavia_gmf(fn) :: LightMesh

    selv, vertex  = GetMesh_GMF(fn)
    vs = Vertex.(mapslices(NTuple{3, Float32}, vertex; dims=2)[:, 1])
    faces = mapslices(NTuple{3, FaceIndex}, selv, dims=2)[:, 1]
     
    return LightMesh(vs, faces)

end

function loadmesh_Lavia_curvo(fn) ::CurvedMesh

    q = GetMeshQuad(fn)
    selvMat = q[1][:, 1:3]
    selvs   = mapslices(NTuple{3}, selvMat; dims=2)[:, 1]
    
    vertexMat = q[2]
    vertices = Vertex.(mapslices(NTuple{3}, vertexMat; dims=2)[:, 1])
    mesh = LightMesh(vertices, selvs) |> cleanupVertices
    
    
    controlMat = q[1][:, 4:6]
    ctrlIndices = mapslices(NTuple{3}, controlMat; dims=2)[:, 1]
    ctrlPts = map.(x -> vertices[x], ctrlIndices)
    
    return CurvedMesh{LightMesh, CFace}(mesh, ctrlPts)
	
end

