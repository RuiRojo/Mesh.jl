#Mesh

```@contents
Pages = ["index.md"]
Depth = 3
```


Lalalal

```@index
Pages = ["index.md"]
Modules = [Mesh]
Order   = [:function, :type]
```

Lololo

Julia block

```julia
x = 23
f(x) :: Int
```

Paquetín para manejar meshes the DAS way.

```@docs
AbstractFace
AbstractMesh
```

# Algunas otras pavadas de prueba

## Autodocs


## Example

```@example

a = 1
b = 2

a + b
```


# Math

Here's some inline maths: ``\sqrt[n]{1 + x + x^2 + \ldots}``.


Here's an equation:

```math
\frac{n!}{k!(n - k)!} = \binom{n}{k}
```

This is the binomial coefficient.

# Autodocs

```@autodocs
Modules = [Mesh]
Order   = [:function, :type]
Private = false
```

Let's recap.


A `Vertex`, holds 3 single precision numbers. 

```julia
v1 = Vertex(0, 0, 0)
v2 = Vertex(1, 0, 0)
v3 = Vertex(0, 1, 0)
```

A Face is build from 3 vertices

```julia
f = Face(v1, v2, v3)
```

The function `vertices` returns the vertices of a face or mesh, and `normal` returns the
normal of a face or triplet of vertices.

An `AbstractFace` requires 

Requiere:

* `vertices`

Implementa:
* `normal`
A 