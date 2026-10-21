# From the repository root: julia --project reproduce.jl
# --skip-shrinking skips ONLY a previously verified coverage certificate.
# --integer selects the independent (slower) integer reconstruction backend.
include("Measurements.jl")
using .SingletPOVMs, BellPolytopes, Serialization, MultiFloats, Printf, SHA, TOML

"""
    reproduce(; verify_shrinking=true, backend=:blas, T=Float64x2,
                storage_path=..., write_report=true, output_dir=@__DIR__)

Verify both measurement-cover certificates and the saved finite local model.
The generated result.toml is an output, never a proof input.
The final visibility retains Alice's exact inradius as a radical expression.
The returned certificate adds visibility_exact; visibility_lower remains the
rational fallback returned by active_set_certificate.
All required geometry, simulation witnesses, and saved atoms are read relative
to this script; output_dir may point to a separate writable report directory.
Use verify_shrinking=false only after separately checking certify_shrinking().
"""
function reproduce(; verify_shrinking = true, backend = :blas, T = Float64x2, storage_path = joinpath(@__DIR__, "s2306_121_v6900.dat"), write_report = true, output_dir = @__DIR__)
    verify_shrinking && certify_shrinking()
    A = measurement_set(; T)
    B = template_set(; T)
    # Certification needs q and exact enclosures, not dense p/o tensors.
    q = SingletPOVMs.orbit_tensor(A, B, false, SingletPOVMs.orbit_dimension(A, B))
    @assert maximum(q) == 32712
    data = certified_povm_data(A, B)
    storage = deserialize(storage_path)
    cert = active_set_certificate(storage, data; v0 = 69 // 100, q, backend)
    exact = exact_alice_visibility(cert)
    println("Certified visibility (exact): ", exact.expression)
    @printf("Decimal display: %.15f\n", Float64((exact.lower + exact.upper) / 2))
    println("Rational fallback: ", cert.visibility_lower)
    println("Face radius: ", cert.radius)
    println("Residual norm upper bound: ", cert.residual_norm_upper)
    println("Analyticity factor: ", cert.nu)
    if write_report
        fraction(x) = "$(numerator(x))/$(denominator(x))"
        report = Dict(
            "format_version" => 2,
            "visibility" => exact.expression,
            "v0" => fraction(cert.v0),
            "visibility_squared_coefficients" => fraction.(collect(exact.squared.c)),
            "field_basis" => ["1", "sqrt(5)", "sqrt((15-3*sqrt(5))/10)", "sqrt(5)*sqrt((15-3*sqrt(5))/10)"],
            "visibility_prefactor" => fraction(exact.prefactor),
            "visibility_enclosure_lower" => fraction(exact.lower),
            "visibility_enclosure_upper" => fraction(exact.upper),
            "visibility_rational_fallback" => fraction(cert.visibility_lower),
            "eta_A" => SingletPOVMs.ALICE_SHRINKING_EXPRESSION,
            "eta_A_squared_coefficients" => fraction.(collect(exact.alice_squared.c)),
            "eta_A_rational_fallback" => fraction(cert.eta[1]),
            "eta_B" => fraction(cert.eta[2]),
            "face_radius" => fraction(cert.radius),
            "residual_upper" => fraction(cert.residual_norm_upper),
            "nu" => fraction(cert.nu),
            "atoms" => length(storage.weights),
            "orbits" => Int(maximum(q)),
            "weight_bits" => cert.weight_bits,
            "coordinate_bits" => data.bits,
            "residual_bits" => cert.residual_bits,
            "radius_bits" => cert.radius_bits,
            "backend" => String(backend),
            "shrinking_verified_this_run" => verify_shrinking,
            "active_set_sha256" => bytes2hex(sha256(read(storage_path))),
            "witness_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "centers_2_witnesses.tsv")))),
            "cover_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "outer_cover.toml")))),
            "targets_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "centers_2_targets.tsv"))))
        )
        mkpath(output_dir)
        open(joinpath(output_dir, "result.toml"), "w") do io
            TOML.print(io, report)
        end
    end
    merge(cert, (visibility_exact = exact,))
end

if abspath(PROGRAM_FILE) == @__FILE__
    all(a -> a in ("--skip-shrinking", "--integer"), ARGS) || error("unknown command-line argument")
    reproduce(; verify_shrinking = !("--skip-shrinking" in ARGS), backend = "--integer" in ARGS ? :integer : :blas)
end
