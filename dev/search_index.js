var documenterSearchIndex = {"docs":
[{"location":"api_reference/#PositiveIntegrators.jl-API","page":"API reference","title":"PositiveIntegrators.jl API","text":"","category":"section"},{"location":"api_reference/","page":"API reference","title":"API reference","text":"CurrentModule = RKOpt","category":"page"},{"location":"api_reference/","page":"API reference","title":"API reference","text":"Modules = [RKOpt]","category":"page"},{"location":"api_reference/#RKOpt.optimize_stability_polynomial-Tuple{Any, Any, Any}","page":"API reference","title":"RKOpt.optimize_stability_polynomial","text":"optimize_stability_polynomial(accuracy_order,\n                              number_of_stages,\n                              spectrum;\n                              dt_min = 0.0,\n                              dt_max = 2.01 * number_of_stages^2 *\n                                       maximum(abs, spectrum),\n                              tol_bisect = 1.0e-9,\n                              tol_feasible = 1.0e-9,\n                              maxiters = 1000,\n                              optimizer = ECOS.Optimizer,\n                              silent = true)\n\nOptimize the stability polynomial of an explicit Runge-Kutta method with a given accuracy_order and number_of_stages for a given spectrum of eigenvalues (a vector of real or complex numbers) using the algorithm of Ketcheson and Ahmadia (2012).\n\nOptional keyword arguments\n\ndt_min = 0.0: Minimum time step size for the bisection.\ndt_max = 2.01 * number_of_stages^2 * maximum(abs, spectrum): Maximum time step size for the bisection.\ntol_bisect = 1.0e-9: Relative tolerance used to terminate the bisection.\ntol_feasible = 1.0e-9: Relative tolerance used to determine the feasibility of the optimization. The problem is considered feasible if the maximal absolute value of the stability polynomial is less than or equal to 1 + tol_feasible at all eigenvalues in the spectrum.\nmaxiters = 1000: Maximum number of iterations for the bisection.\noptimizer = ECOS.Optimizer: Convex optimization solver to be used. You can choose a solver supporting SOCP from the list of supported solvers in the JuMP documentation.\nsilent = true: Whether to suppress the output of the convex optimization solver.\n\nReferences\n\nKetcheson and Ahmadia (2012) Optimal stability polynomials for numerical integration of initial value problems DOI: 10.2140/camcos.2012.7.247\n\n\n\n\n\n","category":"method"},{"location":"code_of_conduct/","page":"Code of conduct","title":"Code of conduct","text":"EditURL = \"https://github.com/NumericalMathematics/RKOpt.jl/blob/main/CODE_OF_CONDUCT.md\"","category":"page"},{"location":"code_of_conduct/#code-of-conduct","page":"Code of conduct","title":"Code of Conduct","text":"","category":"section"},{"location":"code_of_conduct/","page":"Code of conduct","title":"Code of conduct","text":"Contributor Covenant Code of ConductOur PledgeWe as members, contributors, and leaders pledge to make participation in our community a harassment-free experience for everyone, regardless of age, body size, visible or invisible disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education, socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.We pledge to act and interact in ways that contribute to an open, welcoming, diverse, inclusive, and healthy community.Our StandardsExamples of behavior that contributes to a positive environment for our community include:Demonstrating empathy and kindness toward other people\nBeing respectful of differing opinions, viewpoints, and experiences\nGiving and gracefully accepting constructive feedback\nAccepting responsibility and apologizing to those affected by our mistakes, and learning from the experience\nFocusing on what is best not just for us as individuals, but for the overall communityExamples of unacceptable behavior include:The use of sexualized language or imagery, and sexual attention or advances of any kind\nTrolling, insulting or derogatory comments, and personal or political attacks\nPublic or private harassment\nPublishing others' private information, such as a physical or email address, without their explicit permission\nOther conduct which could reasonably be considered inappropriate in a professional settingEnforcement ResponsibilitiesCommunity leaders are responsible for clarifying and enforcing our standards of acceptable behavior and will take appropriate and fair corrective action in response to any behavior that they deem inappropriate, threatening, offensive, or harmful.Community leaders have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, and will communicate reasons for moderation decisions when appropriate.ScopeThis Code of Conduct applies within all community spaces, and also applies when an individual is officially representing the community in public spaces. Examples of representing our community include using an official e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event.EnforcementInstances of abusive, harassing, or otherwise unacceptable behavior may be reported to Hendrik Ranocha. All complaints will be reviewed and investigated promptly and fairly.All community leaders are obligated to respect the privacy and security of the reporter of any incident.Enforcement GuidelinesCommunity leaders will follow these Community Impact Guidelines in determining the consequences for any action they deem in violation of this Code of Conduct:1. CorrectionCommunity Impact: Use of inappropriate language or other behavior deemed unprofessional or unwelcome in the community.Consequence: A private, written warning from community leaders, providing clarity around the nature of the violation and an explanation of why the behavior was inappropriate. A public apology may be requested.2. WarningCommunity Impact: A violation through a single incident or series of actions.Consequence: A warning with consequences for continued behavior. No interaction with the people involved, including unsolicited interaction with those enforcing the Code of Conduct, for a specified period of time. This includes avoiding interactions in community spaces as well as external channels like social media. Violating these terms may lead to a temporary or permanent ban.3. Temporary BanCommunity Impact: A serious violation of community standards, including sustained inappropriate behavior.Consequence: A temporary ban from any sort of interaction or public communication with the community for a specified period of time. No public or private interaction with the people involved, including unsolicited interaction with those enforcing the Code of Conduct, is allowed during this period. Violating these terms may lead to a permanent ban.4. Permanent BanCommunity Impact: Demonstrating a pattern of violation of community standards, including sustained inappropriate behavior,  harassment of an individual, or aggression toward or disparagement of classes of individuals.Consequence: A permanent ban from any sort of public interaction within the community.AttributionThis Code of Conduct is adapted from the [Contributor Covenant][homepage], version 2.0, available at https://www.contributor-covenant.org/version/2/0/codeofconduct.html.Community Impact Guidelines were inspired by Mozilla's code of conduct enforcement ladder.[homepage]: https://www.contributor-covenant.orgFor answers to common questions about this code of conduct, see the FAQ at https://www.contributor-covenant.org/faq. Translations are available at https://www.contributor-covenant.org/translations.","category":"page"},{"location":"license/","page":"License","title":"License","text":"EditURL = \"https://github.com/NumericalMathematics/RKOpt.jl/blob/main/LICENSE\"","category":"page"},{"location":"license/#License","page":"License","title":"License","text":"","category":"section"},{"location":"license/","page":"License","title":"License","text":"MIT LicenseCopyright (c) 2025 Hendrik RanochaPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"EditURL = \"https://github.com/NumericalMathematics/RKOpt.jl/blob/main/CONTRIBUTING.md\"","category":"page"},{"location":"contributing/#Contributing","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"ContributingRKOpt.jl is an open-source project and we are very happy to accept contributions from the community. Please feel free to open issues or submit patches (preferably as pull requests) any time. For planned larger contributions, it is often beneficial to get in contact first, for example via issues.RKOpt.jl and its contributions are licensed under the MIT license (see License). As a contributor, you certify that all your contributions are in conformance with the Developer Certificate of Origin (Version 1.1), which is reproduced below.Developer Certificate of Origin (Version 1.1)The following text was taken from https://developercertificate.org:Developer Certificate of Origin\nVersion 1.1\n\nCopyright (C) 2004, 2006 The Linux Foundation and its contributors.\n1 Letterman Drive\nSuite D4700\nSan Francisco, CA, 94129\n\nEveryone is permitted to copy and distribute verbatim copies of this\nlicense document, but changing it is not allowed.\n\n\nDeveloper's Certificate of Origin 1.1\n\nBy making a contribution to this project, I certify that:\n\n(a) The contribution was created in whole or in part by me and I\n    have the right to submit it under the open source license\n    indicated in the file; or\n\n(b) The contribution is based upon previous work that, to the best\n    of my knowledge, is covered under an appropriate open source\n    license and I have the right under that license to submit that\n    work with modifications, whether created in whole or in part\n    by me, under the same open source license (unless I am\n    permitted to submit under a different license), as indicated\n    in the file; or\n\n(c) The contribution was provided directly to me by some other\n    person who certified (a), (b) or (c) and I have not modified\n    it.\n\n(d) I understand and agree that this project and the contribution\n    are public and that a record of the contribution (including all\n    personal information I submit with it, including my sign-off) is\n    maintained indefinitely and may be redistributed consistent with\n    this project or the open source license(s) involved.","category":"page"},{"location":"#RKOpt.jl","page":"Home","title":"RKOpt.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Julia package RKOpt.jl is work in progress. A mature version of much more functionality is available in the MATLAB-based package RK-Opt.","category":"page"},{"location":"#Example:-Optimize-a-stability-polynomial-for-imaginary-axis-inclusion","page":"Home","title":"Example: Optimize a stability polynomial for imaginary axis inclusion","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"First, we optimize the stability polynomial for imaginary axis inclusion using optimize_stability_polynomial.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using RKOpt: optimize_stability_polynomial\n\nspectrum = im * range(0.0, 1.0, length = 101)\naccuracy_order = 1\nnumber_of_stages = 7\ndt, coefficients = optimize_stability_polynomial(accuracy_order, number_of_stages, spectrum)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Next, we plot the corresponding stability region.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GLMakie\nusing LaTeXStrings\n\nn = 1000\nx = range(-6, 1, length = n)\ny = range(-7, 7, length = 2 * n)\nz = @. x + im * y'\nstab = map(z -> abs(evalpoly(z, coefficients)) <= 1, z)\nfig = Figure()\nax = Axis(fig[1, 1]; aspect = DataAspect(),\n          xlabel = L\"\\mathrm{Re}(z)\", ylabel = L\"\\mathrm{Im}(z)\")\nheatmap!(ax, x, y, stab, colormap = :binary, colorrange = (0.0, 2.0))\nlines!(ax, [extrema(x)...], [0, 0], color = :black, linestyle = :dash)\nlines!(ax, [0, 0], [extrema(y)...], color = :black, linestyle = :dash)\nfig","category":"page"}]
}
