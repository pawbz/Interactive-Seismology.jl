### A Pluto.jl notebook ###
# v0.19.29

#> [frontmatter]
#> title = "Seismology: An Interdisciplinary Science"
#> description = "Some notes on the historical review and introductory material"

using Markdown
using InteractiveUtils

# ╔═╡ 7e939816-34f2-11ee-064e-e3d99343806f
using PlutoUI, PlutoTeachingTools

# ╔═╡ ddecd7e5-4fa5-47a4-a794-acff651f0736
ChooseDisplayMode()

# ╔═╡ 76fda2d0-8178-40d5-bd4c-4a34677f0344
TableOfContents()

# ╔═╡ 08ff2c07-5fe0-41b5-8586-dda51f4c675a
md"""
# Seismology
"""

# ╔═╡ 8002bbb4-3b48-4ab6-a1dc-abd00e4d17a1
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
**the study of mechanical vibrations of the Earth**
	
The science of **seismology** has double feature, it aims simultaneously to obtain
* the infrastructure of the Earth's interior with the aid of **seismic wave phenomena** 
* nature of **earthquake sources** with an ultimate goal of **mitigating** them
    """]))

# ╔═╡ ae4d8b93-0c2c-495a-8c87-fc77521eb1ae
md"# Seismology: Unique Characteristics"

# ╔═╡ 17c11617-d69f-4e39-b58f-1a3426beaabc
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
1) concerned only with the **mechanical properties** and **dynamics** of the Earth
2) **seismic waves** allow investigation of Earth's interior out to the greatest depths than any other branch of geophysics
    * they have **least distortion** compared to any other wave that can be observed after propogating through the Earth
3) contributes to our knowledge of only the **present state** of the Earth's interior
    """]))


# ╔═╡ 8624c99e-88dd-4bce-bf73-c7633f24a2d7
md"# Seismology: Interdisciplinary Science"

# ╔═╡ bc0dd31e-aa00-4fca-8dab-2dfbd4820f0f
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
* medicine: knowledge of anatomy, chemistry of drugs, physics of lasers, mathematics of tomography
* seismology: **geology**, **engineering**, **physics**, **mathematics**, and **computation**
	
    """]))

# ╔═╡ 8c90f08a-f981-4423-bcd7-998a6bc05c3b
Markdown.MD(Markdown.Admonition("note", "",
    [md"""
* most of the theory needed to interpret **seismograms** through efforts of **physicists** and **mathematicians** prior to 1922
* its history is inseparable from the history of great achievements in **continuum mechanics**, **applied mathematics**, and **general wave theory**
    """]))

# ╔═╡ d4d80c02-8d7a-4f38-95f9-2d4f7f5e950b
md"# Seismology"

# ╔═╡ 0d947330-be54-4432-ad05-00838a177682
Markdown.MD(Markdown.Admonition("warning", "Principal Components",
    [md"""
* seismometry and experimental seismology
* theory of seismic fields in the earth (rays, waves, modes, beams)
* seismic sources
* Earth's internal structure 
    """]))

# ╔═╡ 7bd22a44-e063-402e-a107-170d8f866596
md"# Scales in Seismology"

# ╔═╡ 9a07bc8a-02dd-48eb-86df-934bf4ac13c6
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
* **seismic sources**, a factor of $10^{18}$ in equivalent seismic moment
    * smallest detectable micro earthquake to 1960 Chile earthquake
* **seismograph networks**, a factor of $10^6$ in linear dimension
    * 10 meter in engineering surveys to $10000$ km global array observatories
* **ground displacement**, a factor of $10^{11}$
    * smallest detectable ground displacement of $10^{-10}$ m (compare with diameter of the hydrogen atom) to ≈$10$ m of slip on a major fault during the earthquake    """]))

# ╔═╡ 3dddccb3-3559-49d9-8b83-243501fb3732
md"""
# Seesaws in Seismology
"""

# ╔═╡ 00d8ac09-7674-40f3-9cd8-0314102e3e2a
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
effects of sources and medium are strongly **coupled** in seismology
* our knowledge of the **seismic sources** and the **Earth medium** have advanced in a see-saw fashion, i.e., 
* if sources are better understood at a given stage, they are constrained to better understand the medium, and then the improved knowledge of the medium is used to revise our knowledge of the source
    """]))

# ╔═╡ cd38117c-9eb4-4342-8444-18ad400961cf
md"# Seismology"

# ╔═╡ 010c21f4-3f60-4f9a-b2d3-6ab6616aca83
Markdown.MD(Markdown.Admonition("warning", "Stimuli of Rapid Growth",
    [md"""
* occurrence of major devastating earthquakes
* exploration (oil and gas, minerals)
* advances in mathematics and theoretical physics
* breakthroughs in sensing or computation technology
* planetary exploration
    """]))

# ╔═╡ dcc341d3-8d91-43f4-9848-d919f10518cb
md"# Seismology and Big Data"

# ╔═╡ e30f6c5c-e016-4a23-93f3-0f535c6082fe
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
exponential growth in better quality data
* 1901: 25 seismic stations
* 1940: 250 seismic stations around the globe
* 2023: 10000 seismic stations around the globe
    """]))

# ╔═╡ f142b945-7baf-4737-a14e-24fd917ea47b
Markdown.MD(Markdown.Admonition("warning", "Current Holistic Statge",
    [md"""
* imaging the Earth's interior and seismic sources will all available data
    """]))

# ╔═╡ ac9d894c-d83a-4005-b2c7-d7e3f42f9b6b
md"""
# Seismologists are Handicapped
"""

# ╔═╡ a65d7d91-f18b-4202-b48a-1c674bf3b6a2
Markdown.MD(Markdown.Admonition("danger", "",
    [md"""
* meteorologists can put their sensors in the eye of a hurricane, we cannot put a sensor in the focal region of an earthquake
* heavy dependency on the far-field measurements
* **physical theory** as to what happens at the earthquake source is missing
    """]))

# ╔═╡ 30c7c3c2-9f04-46b8-aff7-5d5420f70279
Markdown.MD(Markdown.Admonition("danger", "",
    [md"""
* most of the current graduate students engaged in _computer simulation_ games and/or _seismic-data processing_ challenges, with insufficient efforts to construct **new theoretical** physical models
* we need endeavors to come up with **new mathematical weapons** to understand the **nonlinear dynamics**
    """]))

# ╔═╡ 0273dc65-f4ed-46b3-884d-66cf9e4e3401
md"""
# Some Technical Remarks
"""

# ╔═╡ a4c39dca-0d55-4d7a-9715-7e06d2c52cc6
Markdown.MD(Markdown.Admonition("danger", "",
    [md"""
	
A more technical definition:  
* mostly a study where the spatial fluctuations in particle displacements, strains, and stresses have wavelengths much larger than the amplitudes of the particle displacements in the Earth's medium
* in other words, we don't use the equation of motion in its *strict form* (Aki and Richards, Ch.2)
    """]))

# ╔═╡ bbab7f7a-cf8e-47f4-8484-339bc732b3ee
md"# Historical Review"

# ╔═╡ 8150e6db-4619-4dfa-8a94-8b39a1299a6a
md"""
## Seismometry
"""

# ╔═╡ 6fb5676c-7cf2-4c41-a2b5-95963b5e3ce9
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
    a seismograph for an earth scientist what the telescope is to the astronomer; a tool for peering into inaccessible regions -- Ari Ben-Menahem
    """]))

# ╔═╡ 91afe5f3-af6a-4aeb-8522-9fddb33b0ccb
md"""
### 1819
Observations of faulting associated with Kutch earthquake in India
"""

# ╔═╡ b43748a6-ddbb-4b65-a3e3-3a090344c9c6
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
major earthquakes continue to serve as milestones on the road of progress in seismology
    """]))

# ╔═╡ bebd002d-6d6b-4205-9aa1-77a69bb7b75e
md"""
### 1841 James David Forbes
* the first seismometer 
"""

# ╔═╡ 3ef6f004-cd6d-450b-858a-755dafc1e7b0
md"""
### 1857 Robert Mallet
* an engineer who laid the foundation of instrumental seismology
* first systematic attempt to apply physical principles to earthquake effects
* first world seismicity map
"""

# ╔═╡ f45f25cd-91e1-4f38-aebd-b354a421b0f9
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
Mallet, the first true seismologist, used bowls of mercury at varying distances from sources to measure the propogation of waves!
    """]))

# ╔═╡ 428d6419-cf4f-4a66-b189-84e94887fe44
md"""
### 1884 Rossi-Forel, J. F. J Schmidt
* a scale for earthquake effects
"""

# ╔═╡ e88b8fc6-9add-42de-9e57-ab5c9b96bdfb
md"""
### 1889 Ernst von Rebeur-Paschwitz
* first measurement of teleseismic (Japanese) earthquake
"""

# ╔═╡ 279c7894-71a4-44b2-9329-3eb55f6d5e19
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
birth of teleseismic seismograms: earthquakes are no longer a local affair of cities but became a global phenomenon -- this was the **first revolution**
    """]))

# ╔═╡ 4c632977-3fe1-49c5-8d0a-2ea293c69ee1
Markdown.MD(Markdown.Admonition("tems", "",
    [md"""
When first seismogram was measured in early 1880s, seismologists were puzzled why it  lasted long and why there are oscillations after the arrivals or P and S waves!
    """]))

# ╔═╡ f835861f-5dac-408f-b3ce-f7e30cc55bdb
md"""
### 1894 Milne and his associates
* compact seismograph system, and seismology emerged as a quantitative science
"""

# ╔═╡ 3cd8c1ce-fe1b-48e1-a27a-3aa9e0a6a736
md"""
### 1895 F. Omori
* a law for aftershock time series
"""

# ╔═╡ ce64efef-def6-40f2-8795-7c51cf54fed5
md"""
### 1897 R. D. Oldham
* identification of three types of waves (P, S and R) in seismograms
"""

# ╔═╡ 1c12869d-9bae-446b-8189-9019663ff4d3
md"""
### 1900 Emil Wiechert
* three-component mechanical seismograph system
"""

# ╔═╡ 512b93dd-bbb0-4610-95fb-c7c65d2a19e1
md"""
### 1906 
* observed faulting and slip for California earthquake
"""

# ╔═╡ ed2bc3a9-f898-4fe3-b013-8636dd4803d3
md"""
### 1910 Boris Borisovich Golitzin
* first electromagnetic seismograph
"""

# ╔═╡ a1f5916d-edf2-444d-90ef-3ab83bb00e0e
md"""
### 1935 Hugo Benioff
* measure a component of ground _strain_; a strain seismograph measures variation in displacement between two points in space
"""

# ╔═╡ bf68f909-c887-4637-a792-c697666c699f
md"""
### 1942 B. Gutenberg, C. Richter
* first empirical observations between earthquake magnitude, intensity, energy, and frequency of occurrence
"""

# ╔═╡ f5a7a700-c66a-424d-8c8c-d8f4c1e3d17e
md"""
### 1946
* age of nuclear testing; the use of nuclear explosions greatly enhanced the capabilities of seismic studies of the Earth's interior
"""

# ╔═╡ 608259ba-fb7d-41b4-8b6b-33b0d2a0d377
md"""
### 1964--1969
* large seismometer arrays came into vogue
"""

# ╔═╡ 49b41407-bb0d-46e5-a020-1443f8c78b16
md"""
### 1969 Apollo passive seismic experiment
* beginning of planetary seismology
"""

# ╔═╡ d4a1417c-d424-44c4-bcdd-4b601466200b
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
When the first seismogram for the Moon was obtained in 1969, seismologists were again  puzzled by the great length of time for which oscillations continued!
    """]))

# ╔═╡ b2a52562-9bd9-4d8a-91ba-de0f4a2b1007
md"""
### 1978
* ocean-bottom seismographs
"""

# ╔═╡ 7e9bf8ee-1ce4-4833-8575-4d6384e4cc07
md"""
### 1982 Global Digital Seismograph Network
* a worldwide network of more than 100 stations
"""

# ╔═╡ 49934ec7-0e10-4224-a5af-d0630c8c5d6c
md"""
### 2011 Distributed Acoustic Sensing (DAS)
* Rayleigh scattering-based  DAS systems use fiber optic cables to provide distributed strain sensing
"""


# ╔═╡ c0a52c2f-61b6-4693-9d7b-99f43c56225d
md"""
## Theory of Seismic Fields in the Earth
"""

# ╔═╡ 8bc517c9-1657-4329-8cc7-7a06699864d7
md"""
### 1660 Robert Hooke
* stated one-dimensional linear stress-strain relationship, thus laying the foundation for theory of elasticity
"""

# ╔═╡ 6538053f-254c-4d14-b24a-0d9a834700b3
md"""
### 1807 Thomas Young
* recognize shear as elastic strain and defined modulus of elasticity
"""

# ╔═╡ d6ec37b9-d305-47f1-922c-b7c99cede638
md"""
### 1821 to 1831 C. Navier, A. Cauchy, S. D. Poisson
* fundamental equations of linear elastodynamics
* Cauchy's stress-strain relations
* existence of compressional and shear waves in elastic solids
* created confusion in the wave theory of light if two types of light waves should be visible
"""

# ╔═╡ f9829464-b5c5-4279-a478-97a2a732b79e
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
Navier was a disciple of Fourier and a professor of mechanics; he was the first one to derive elastodynamic displacement equations
    """]))

# ╔═╡ e1f9dec3-d4c3-4651-b6ee-bf6dfa1f9d75
md"""
### 1829 to 1850 S. D. Poisson, G. G. Stokes, Lord Kelvin, H. Lamb, G. Lame
* theoretical studies on the vibrations of elastic bodies
* Stokes fundemental equation on viscous fluids
"""

# ╔═╡ 782fc90f-9aef-46d6-ba85-3c88a7d98aa1
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
even though Poisson discovered longitudinal and transverse waves in 1828, they were only first observed in seismograms by Oldham in 1897!
    """]))

# ╔═╡ 048baaad-be20-4b31-807f-1e999bd97a6f
md"""
### 1872 E. Betti
* Reciprocity relations leading later to the representation theorem of elastodynamics
"""

# ╔═╡ 0f96094d-79bb-4bfd-b3a6-acf0b1e5065d
md"""
### 1877 E. B. Christoffel
* plane-wave phase velocities in general anisotropic media
"""

# ╔═╡ bc237ace-106d-4fda-811e-fbd0b22e669c
md"""
### 1885 Lord Rayleigh
* a homogeneous elastic substance can accommodate a third wave at its boundary
* Rayleigh waves are inhomogeneous waves

"""

# ╔═╡ 206e27b3-09b6-4102-918c-9bd82c031824
md"""
### 1898 T. Bromwich
* studied the influence of gravity on elastic waves
"""

# ╔═╡ facc0073-7a09-43bc-a631-3996fc70eacc
md"""
### 1899 C. G. Knott
* refraction and reflection coefficients of plane seismic waves at planar discontinuities 
"""

# ╔═╡ ac9ed2ee-0f82-45f5-b797-70117714ba90
md"""
### 1904 H. Lamb
* first synthetic seismogram: a mathematical model of an earthquake in a half-space configuration, known as Lamb's problem
"""

# ╔═╡ 8b15fcc3-947f-4fbd-b78b-8575d2281c92
md"""
### 1904 A. E. H. Love
* 3-D elastodynamic representation theorem 
"""

# ╔═╡ 807d8372-3bbe-49c6-95da-a27527613ed1
md"""
### 1905 R. Becker
* theory of standard-linear solids with relaxation, as a physical model for anelasticity
"""

# ╔═╡ 442b48e4-feaf-444e-b8ce-02c43995a73b
md"""
### 1910 H. Benndrof, K. Zoppritz
* ray theory for seismic signal propagation in spherically symmetric Earth models
"""

# ╔═╡ 844cc4b3-c0b0-4e42-8d55-219798ebd762
md"""
### 1910 M. R. Frechet
* functional calculus (the "Frechet" derivative)
"""

# ╔═╡ 13e8e5b8-4d6e-4694-b359-dc6bd75be105
md"""
### 1911 A. E. H. Love
* _Love waves_ are transversely polarized waves, which are not included in the theories of Rayleigh and Lamb
"""

# ╔═╡ d80f6fc8-e523-45b2-abf4-cc282f252cf5
md"""
### 1912 H. Lamb
* group and phase velocity of surface waves and the interpretation of seismograms
"""

# ╔═╡ 7d6ae8e4-1e17-4cb4-bb7b-64cbdd7e155f
md"""
### 1920 Harold Jeffreys
* advanced mathematical and statistical methods into seismic data analysis
"""

# ╔═╡ deb743fa-0c2f-4340-a41c-9c7f7757c0b5
md"""
### 1920 L. M. Hoskins
* studied free oscillations of gravitating radially inhomogeneous spheres
"""

# ╔═╡ c277e885-8e1e-4819-ab3a-1beb34344d32
md"""
### 1926 M. Born
* the Born approximation of wavefields
"""

# ╔═╡ 8a7e51d4-8fb1-4a2f-8b88-fa4730c2eab2
md"""
### 1930 N. Wiener
* tools for computerized signal analysis 
"""

# ╔═╡ 7d011d40-fb4a-446a-9cc2-598fb13a7e0f
md"""
### 1949 R. Stoneley
* effect of anisotropy on elastic surface waves
"""

# ╔═╡ 3a6544b2-2094-4868-ba16-a79cc2feebe8
md"""
### 1953 P. M. Morse, H. Feshbach
* elastodynamic integral representation theorem 
"""

# ╔═╡ 3e726b0d-b7ff-46e3-8200-ffb8dd345d3c
md"""
### 1955
* Fourier-transform methods into analysis of surface-wave dispersion and attenuation
"""

# ╔═╡ 216e9eb0-6cfa-4bcb-a422-6a3a593e7064
md"""
### 1960--1975
* age of digital computers
* computer-generated synthetic seismograms
* huge amounts of seismological data can be processed in a relatively short time
* asymptotic wave theories in vertically in homogeneous media
"""

# ╔═╡ 47cdf739-6926-4908-9e4b-7d04228ee2eb
md"""
### 1970--1985
* asymptotic ray theory, ray tracing in complex media, shadow zones
* finite-difference and finite-element methods 
* theory of inversion and resolution of gross Earth data
* imaging in inhomogeneous elastic media
* Kirchoff-Helmholtz integral formulation; Kirchoff migration
"""

# ╔═╡ 2ddc2428-b478-4a78-bf70-2cb2a796ca25
md"""
### 1984--
* back-projection of travel-time delays along ray paths
* seismic-wave tomographic inversion
* diffraction tomography in exploration seismology
* full waveform inversion
"""

# ╔═╡ 9ac9cc91-9c77-4103-9381-4191486de3f5
md"""
## Earth's Internal Structure
"""

# ╔═╡ c378ef72-7a7e-4edc-956a-b9499c389534
md"""
### 1798 Henry Cavendish
* determined mean density of the Earth using Newton's law of universal gravitation
"""

# ╔═╡ f6857c5c-baad-4e3c-833c-5fb7a8ec2022
md"""
### 1863 Lord Kelvin
* used tidal observations to estimate that the mean rigidity of the Earth exceeds that of steel
"""

# ╔═╡ 90ad480f-38bc-4c17-a296-ef2a1a72b2e3
md"""
### 1888 A. Schmidt
* argues in general that wave velocity increases with depth, and raypaths will be curved and convave upward toward the Earth surface
"""

# ╔═╡ d20f70e1-8fbb-43ee-846d-1bd1e36867ac
md"""
### 1897 Emil Wiechert
* numerical details of an Earth model consisting of a core of radius almost 5000 km
* conjecture that the inner core is mettalic
"""

# ╔═╡ 53e44c54-ee0f-4180-96f9-0e9e85b24a94
md""" ### 1910 R. D. Oldham
* presence of central core with a radius of 1600 km after observing substantially delayed P waves 
* inferred that the Earth contains a central region with lower wave velocity than the surrounding shell (the mantle)
* no S waves were observed below the mantle

"""

# ╔═╡ 28f35d0b-f406-4a98-a448-39674262bdc9
md"""
### 1909 K. Zoppritz
* first traveltime inversion: calculation of travel-time tables for seismic phases and analytical solutions for Abel-type integral equations to determine wave velocity as a function of depth
"""

# ╔═╡ 45059635-edaf-40b7-8e89-b97551b9ec95
md"""
### 1909 Andrija Mohorovicic
* evidence of a sharp increase in the P-wave velocity 54 km below the Earth's surface
"""

# ╔═╡ 5b2666f9-a442-447c-8265-b3df07b639d4
md"""
### 1914 Beno Gutenberg
* found a low-velocity zone at a depth of 2900 km
"""

# ╔═╡ d7c5a81c-0f82-4361-b25e-2eb4deb30e9f
md"""
### 1921 E. Meissner
* observed dispersion of surface waves in the Earth crust
"""

# ╔═╡ cd4f205b-2bec-450e-be48-72508ce04f94
md"""
### 1923 E. D. Williamson and L. H. Adams
* use seismological data, values of Earth's mass and moment of inertia to estimate density gradients of the Earth
"""

# ╔═╡ 97852f31-2e15-464d-9c9f-fb5fb2a2db52
md"""
### 1930 V. Conard H. Jeffreys
* three-layer crustal models using near-field headwaves ($P_g$, $S_g$, $P^*$, $S^*$)
"""

# ╔═╡ bbf6ac3b-3280-47fe-a490-e4c0e42a7f1d
md"""
### 1936 Inge Lehmann
* existence of an inner core by observing amplitudes of P waves between angular distances of $105°$ to $142°$
"""

# ╔═╡ 95b2372f-b042-4191-89c8-7da386502afa
md"""
### 1939 H. Jeffreys
* apply Airy theory of diffraction near a caustic to understand wave diffraction by the Earth's core
"""

# ╔═╡ 09c9ce19-e085-4422-b8f8-8034b0dd18d4
md"""
### 1942 K. E. Bullen
* classification of Earth's interior into a number of shells
"""

# ╔═╡ 2063c6e1-11fe-4f61-96cd-89f6277e7b8f
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
computers enabled sophisticated forward and inverse computational schemes; this was  the second revolution 1950--1955
    """]))

# ╔═╡ ccb23061-b4d7-4e62-b1f3-e4d5f6c845d7
md"""
### 1969--1977
* rays in anisotropic media; evidence of anisotropy in crust and mantle
"""

# ╔═╡ ec5d6c5e-c61e-468e-bb9a-366409a5fb3b
md"""
### 1980 --
* ever-increasing quantity and quality of seismic data, deviations from isotropy, spherical symmetry, and pure elasticity were observed with more precision

* 3-D images of seismic velocities and anisotropy with detailed velocity structure near major internal boundaries
"""

# ╔═╡ 4dce2716-c51d-4cc2-892e-94d37b1077fc
md"""
## Seismic Sources
"""

# ╔═╡ 67de6071-c446-4740-8e6c-b137d5350443
md"""
### 1761 John Michell
* explosive theory of earthquakes
* earthquakes originate within the Earth, and wave spread out from the source throughout the Earth's interior
"""

# ╔═╡ d5a3f768-7bcb-4d89-a3b2-e75041faf393
md"""
### 1906 H. F. Reid
* studied geodetic measurements before and after an earthquakee rupture, and developed elastic rebound theory
"""

# ╔═╡ fec411a1-e927-4f83-8426-9f2fd4836f4a
md"""
### 1907 Vito Volterra
* integral representation of seismic sources using Betti's reciprocity theorem
"""

# ╔═╡ cd8ec874-e24d-449f-9549-6cd5ba32c461
md"""
### 1917 J. Shida
* observed regularities in the distribution of polarities of the initial P-wave motion
* epicentral region is divided into four parts by two perpendicular lines intersecting at the epicenter
"""

# ╔═╡ b1535bcf-1023-4f94-bec5-12683479d321
md"""
### 1922 H. H. Turner
* existance of deep earthquakes
"""

# ╔═╡ 1f0dc467-3e93-4baf-a1cb-0912ba0dae82
md"""
### 1923 Hiroshi Nakano
* double-couple model: observed patterns of initial motions are explainable using classical Stokes-Love solution
"""

# ╔═╡ 64b7e782-827e-4caf-84ea-a8a6fcf93201
md"""
### 1926 H. Jeffreys
* first use of wave amplitudes from seismograms to estimate earthquake energy release
"""

# ╔═╡ 52c0be56-b35c-4b29-b6ca-99a2c7ec7402
md"""
### 1952 B. Gutenberg, C. Richter
* surface-wave magnitude: maxiumum amplitude of the crustal surface waves having a specific period
"""

# ╔═╡ 53c70c64-4a9e-411d-a6e4-5a3bd793e387
md"""
### 1958 J. A. Steketee
* equivalence theorem, stating that the displacement field produced by dislocation sources in an elastic body equals that produced by double couple forces
"""

# ╔═╡ 0aebcf33-15c2-4a00-a13c-91e582baf11d
md"""
### 1960
* Chile earthquake was studied comprehensively
"""

# ╔═╡ b36647c9-248d-4047-840b-531401ca3c70
md"""
### 1960--1980
* propagating rupture over a causative fault is revealed
* kinematic dislocation model 
* energy budget for major earthquakes
"""

# ╔═╡ 598d003d-0932-4fcc-988c-ecd69537899c
md"""
### 1989 --
* detailed study of the geometry of the rupture using diversified set of data: near-field and far-field seismograms, GPS measurements, surface ruptures, aftershock spatial distribution, and focal mechanism

* results obtained so far point to complexity of seismic sources
"""

# ╔═╡ 710f75b9-c581-4134-82db-011644c5b81e
md"""
## Resources
[^Paper]: Ari Ben-Menahem; A concise history of mainstream seismology: Origins, legacy, and perspectives. Bulletin of the Seismological Society of America
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoTeachingTools = "~0.2.12"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "2af2aaef1ea67811bc114f2d0bf7f61179b912c7"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "d730914ef30a06732bdd9f763f6cc32e92ffbff1"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "6a125e6a4cb391e0b9adbd1afa9e771c2179f8ef"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.23"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "8c57307b5d9bb3be1ff2da469063628631d4d51e"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.21"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    DiffEqBiologicalExt = "DiffEqBiological"
    ParameterizedFunctionsExt = "DiffEqBase"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    DiffEqBase = "2b5f629d-d688-5b77-993f-72d75c75574e"
    DiffEqBiological = "eb300fae-53e8-50a0-950c-e21f52c2b7e0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "4b2e829ee66d4218e0cef22c0a64ee37cf258c29"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "45f9e1b6f62a006a585885f5eb13fc22554a8865"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.12"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "1e597b93700fa4045d7189afa7c004e0584ea548"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═ddecd7e5-4fa5-47a4-a794-acff651f0736
# ╠═76fda2d0-8178-40d5-bd4c-4a34677f0344
# ╟─08ff2c07-5fe0-41b5-8586-dda51f4c675a
# ╟─8002bbb4-3b48-4ab6-a1dc-abd00e4d17a1
# ╟─ae4d8b93-0c2c-495a-8c87-fc77521eb1ae
# ╟─17c11617-d69f-4e39-b58f-1a3426beaabc
# ╟─8624c99e-88dd-4bce-bf73-c7633f24a2d7
# ╟─bc0dd31e-aa00-4fca-8dab-2dfbd4820f0f
# ╟─8c90f08a-f981-4423-bcd7-998a6bc05c3b
# ╟─d4d80c02-8d7a-4f38-95f9-2d4f7f5e950b
# ╟─0d947330-be54-4432-ad05-00838a177682
# ╟─7bd22a44-e063-402e-a107-170d8f866596
# ╟─9a07bc8a-02dd-48eb-86df-934bf4ac13c6
# ╟─3dddccb3-3559-49d9-8b83-243501fb3732
# ╟─00d8ac09-7674-40f3-9cd8-0314102e3e2a
# ╟─cd38117c-9eb4-4342-8444-18ad400961cf
# ╟─010c21f4-3f60-4f9a-b2d3-6ab6616aca83
# ╟─dcc341d3-8d91-43f4-9848-d919f10518cb
# ╟─e30f6c5c-e016-4a23-93f3-0f535c6082fe
# ╟─f142b945-7baf-4737-a14e-24fd917ea47b
# ╟─ac9d894c-d83a-4005-b2c7-d7e3f42f9b6b
# ╟─a65d7d91-f18b-4202-b48a-1c674bf3b6a2
# ╟─30c7c3c2-9f04-46b8-aff7-5d5420f70279
# ╟─0273dc65-f4ed-46b3-884d-66cf9e4e3401
# ╟─a4c39dca-0d55-4d7a-9715-7e06d2c52cc6
# ╟─bbab7f7a-cf8e-47f4-8484-339bc732b3ee
# ╟─8150e6db-4619-4dfa-8a94-8b39a1299a6a
# ╟─6fb5676c-7cf2-4c41-a2b5-95963b5e3ce9
# ╟─91afe5f3-af6a-4aeb-8522-9fddb33b0ccb
# ╟─b43748a6-ddbb-4b65-a3e3-3a090344c9c6
# ╟─bebd002d-6d6b-4205-9aa1-77a69bb7b75e
# ╟─3ef6f004-cd6d-450b-858a-755dafc1e7b0
# ╟─f45f25cd-91e1-4f38-aebd-b354a421b0f9
# ╟─428d6419-cf4f-4a66-b189-84e94887fe44
# ╟─e88b8fc6-9add-42de-9e57-ab5c9b96bdfb
# ╟─279c7894-71a4-44b2-9329-3eb55f6d5e19
# ╟─4c632977-3fe1-49c5-8d0a-2ea293c69ee1
# ╟─f835861f-5dac-408f-b3ce-f7e30cc55bdb
# ╟─3cd8c1ce-fe1b-48e1-a27a-3aa9e0a6a736
# ╟─ce64efef-def6-40f2-8795-7c51cf54fed5
# ╟─1c12869d-9bae-446b-8189-9019663ff4d3
# ╟─512b93dd-bbb0-4610-95fb-c7c65d2a19e1
# ╟─ed2bc3a9-f898-4fe3-b013-8636dd4803d3
# ╟─a1f5916d-edf2-444d-90ef-3ab83bb00e0e
# ╟─bf68f909-c887-4637-a792-c697666c699f
# ╟─f5a7a700-c66a-424d-8c8c-d8f4c1e3d17e
# ╟─608259ba-fb7d-41b4-8b6b-33b0d2a0d377
# ╟─49b41407-bb0d-46e5-a020-1443f8c78b16
# ╟─d4a1417c-d424-44c4-bcdd-4b601466200b
# ╟─b2a52562-9bd9-4d8a-91ba-de0f4a2b1007
# ╟─7e9bf8ee-1ce4-4833-8575-4d6384e4cc07
# ╟─49934ec7-0e10-4224-a5af-d0630c8c5d6c
# ╟─c0a52c2f-61b6-4693-9d7b-99f43c56225d
# ╟─8bc517c9-1657-4329-8cc7-7a06699864d7
# ╟─6538053f-254c-4d14-b24a-0d9a834700b3
# ╟─d6ec37b9-d305-47f1-922c-b7c99cede638
# ╟─f9829464-b5c5-4279-a478-97a2a732b79e
# ╟─e1f9dec3-d4c3-4651-b6ee-bf6dfa1f9d75
# ╟─782fc90f-9aef-46d6-ba85-3c88a7d98aa1
# ╟─048baaad-be20-4b31-807f-1e999bd97a6f
# ╟─0f96094d-79bb-4bfd-b3a6-acf0b1e5065d
# ╟─bc237ace-106d-4fda-811e-fbd0b22e669c
# ╟─206e27b3-09b6-4102-918c-9bd82c031824
# ╟─facc0073-7a09-43bc-a631-3996fc70eacc
# ╟─ac9ed2ee-0f82-45f5-b797-70117714ba90
# ╟─8b15fcc3-947f-4fbd-b78b-8575d2281c92
# ╟─807d8372-3bbe-49c6-95da-a27527613ed1
# ╟─442b48e4-feaf-444e-b8ce-02c43995a73b
# ╟─844cc4b3-c0b0-4e42-8d55-219798ebd762
# ╟─13e8e5b8-4d6e-4694-b359-dc6bd75be105
# ╟─d80f6fc8-e523-45b2-abf4-cc282f252cf5
# ╟─7d6ae8e4-1e17-4cb4-bb7b-64cbdd7e155f
# ╟─deb743fa-0c2f-4340-a41c-9c7f7757c0b5
# ╟─c277e885-8e1e-4819-ab3a-1beb34344d32
# ╟─8a7e51d4-8fb1-4a2f-8b88-fa4730c2eab2
# ╟─7d011d40-fb4a-446a-9cc2-598fb13a7e0f
# ╟─3a6544b2-2094-4868-ba16-a79cc2feebe8
# ╟─3e726b0d-b7ff-46e3-8200-ffb8dd345d3c
# ╟─216e9eb0-6cfa-4bcb-a422-6a3a593e7064
# ╟─47cdf739-6926-4908-9e4b-7d04228ee2eb
# ╟─2ddc2428-b478-4a78-bf70-2cb2a796ca25
# ╟─9ac9cc91-9c77-4103-9381-4191486de3f5
# ╟─c378ef72-7a7e-4edc-956a-b9499c389534
# ╟─f6857c5c-baad-4e3c-833c-5fb7a8ec2022
# ╟─90ad480f-38bc-4c17-a296-ef2a1a72b2e3
# ╟─d20f70e1-8fbb-43ee-846d-1bd1e36867ac
# ╟─53e44c54-ee0f-4180-96f9-0e9e85b24a94
# ╟─28f35d0b-f406-4a98-a448-39674262bdc9
# ╟─45059635-edaf-40b7-8e89-b97551b9ec95
# ╟─5b2666f9-a442-447c-8265-b3df07b639d4
# ╟─d7c5a81c-0f82-4361-b25e-2eb4deb30e9f
# ╟─cd4f205b-2bec-450e-be48-72508ce04f94
# ╟─97852f31-2e15-464d-9c9f-fb5fb2a2db52
# ╟─bbf6ac3b-3280-47fe-a490-e4c0e42a7f1d
# ╟─95b2372f-b042-4191-89c8-7da386502afa
# ╟─09c9ce19-e085-4422-b8f8-8034b0dd18d4
# ╟─2063c6e1-11fe-4f61-96cd-89f6277e7b8f
# ╟─ccb23061-b4d7-4e62-b1f3-e4d5f6c845d7
# ╟─ec5d6c5e-c61e-468e-bb9a-366409a5fb3b
# ╟─4dce2716-c51d-4cc2-892e-94d37b1077fc
# ╟─67de6071-c446-4740-8e6c-b137d5350443
# ╟─d5a3f768-7bcb-4d89-a3b2-e75041faf393
# ╟─fec411a1-e927-4f83-8426-9f2fd4836f4a
# ╟─cd8ec874-e24d-449f-9549-6cd5ba32c461
# ╟─b1535bcf-1023-4f94-bec5-12683479d321
# ╟─1f0dc467-3e93-4baf-a1cb-0912ba0dae82
# ╟─64b7e782-827e-4caf-84ea-a8a6fcf93201
# ╟─52c0be56-b35c-4b29-b6ca-99a2c7ec7402
# ╟─53c70c64-4a9e-411d-a6e4-5a3bd793e387
# ╟─0aebcf33-15c2-4a00-a13c-91e582baf11d
# ╟─b36647c9-248d-4047-840b-531401ca3c70
# ╟─598d003d-0932-4fcc-988c-ecd69537899c
# ╟─710f75b9-c581-4134-82db-011644c5b81e
# ╠═7e939816-34f2-11ee-064e-e3d99343806f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
