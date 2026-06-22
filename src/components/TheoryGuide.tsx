/**
 * Theoretical Scattering Physics Guide
 * High-quality educational resources detailing Dirac equations,
 * Legendre partial-wave summations, and electrostatic screening.
 */

import React from "react";
import { BookOpen, HelpCircle, Shield, Atom, Sparkles } from "lucide-react";

export default function TheoryGuide() {
  return (
    <div id="theory-guide-section" className="bg-[#12131a] border border-[#222430] rounded-xl p-6 shadow-[0_4px_25px_rgba(0,0,0,0.5)] flex flex-col gap-5 font-sans leading-relaxed text-[#94a3b8]">
      <div className="border-b border-[#222430] pb-3 select-none shrink-0">
        <h2 className="text-base font-bold text-white flex items-center gap-2 font-serif tracking-wider">
          <BookOpen className="w-5 h-5 text-[#c5a059] animate-pulse" />
          Dirac Partial-Wave Atom Elastic Scattering Theory
        </h2>
        <p className="text-xs text-[#94a3b8] mt-1">
          A pedagogical outline of the physical equations resolved inside the ELSEPA simulator engine.
        </p>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-6 text-xs leading-normal">
        {/* Card 1: DCS */}
        <div className="flex flex-col gap-2.5">
          <div className="flex items-center gap-1.5 text-white font-serif font-bold text-[13px] tracking-wide">
            <Atom className="w-4.5 h-4.5 text-[#c5a059]" />
            1. Differential Cross Section (DCS)
          </div>
          <p className="text-[#94a3b8]">
            The differential cross section <code className="bg-[#090a0f] px-1 py-0.5 rounded font-mono text-[#dfba73] border border-[#222430]">dσ/dΩ</code> represents the spatial density distribution of scattered projectiles (electrons or positrons) per unit solid angle. In quantum mechanics, it is mathematically calculated from the complex scattering amplitudes <code className="bg-[#090a0f] px-1 py-0.5 rounded font-mono text-[#dfba73] border border-[#222430]">f(θ)</code> (spin-up) and <code className="bg-[#090a0f] px-1 py-0.5 rounded font-mono text-[#dfba73] border border-[#222430]">g(θ)</code> (spin-flip):
          </p>
          <div className="bg-[#090a0f] p-4 rounded-lg font-mono text-center text-[#c5a059] font-bold border border-[#222430] my-2 text-[14px]">
            dσ/dΩ = |f(θ)|² + |g(θ)|²
          </div>
          <p className="text-[#94a3b8]">
            DCS yields units of <code className="bg-[#090a0f] px-1 py-0.5 rounded font-mono text-[#dfba73] border border-[#222430]">a₀²/sr</code> (Bohr radius squared per steradian) or <code className="bg-[#090a0f] px-1 py-0.5 rounded font-mono text-[#dfba73] border border-[#222430]">cm²/sr</code>.
          </p>
        </div>

        {/* Card 2: Dirac vs Schrodinger */}
        <div className="flex flex-col gap-2.5">
          <div className="flex items-center gap-1.5 text-white font-serif font-bold text-[13px] tracking-wide">
            <Shield className="w-4.5 h-4.5 text-purple-400" />
            2. Dirac Relativistic Formulation
          </div>
          <p className="text-[#94a3b8]">
            While Schrödinger's mechanics is sufficient for light items, heavy elements (high <code className="bg-[#090a0f] px-1 py-0.5 rounded text-[#dfba73] border border-[#222430] font-mono">Z</code> gold, uranium) and high energy projectiles (<code className="bg-[#090a0f] px-1 py-0.5 rounded font-mono text-[#dfba73] border border-[#222430]">E &gt; 10 keV</code>) demand Einsteinian relativity. The projectile electron's speed approaches significant fractions of <code className="bg-[#090a0f] px-1 py-0.5 rounded text-[#dfba73] border border-[#222430] font-mono">c</code>, inducing mass increase and strong spin-orbit coupling.
          </p>
          <p className="text-[#94a3b8]">
            ELSEPA solves the **radial Dirac equation** to produce exact phase shifts ($d_l$ and $e_l$), giving extremely high fidelity angular cross-section shapes.
          </p>
        </div>

        {/* Card 3: Phase shifts */}
        <div className="flex flex-col gap-2.5">
          <div className="flex items-center gap-1.5 text-white font-serif font-bold text-[13px] tracking-wide">
            <Sparkles className="w-4.5 h-4.5 text-amber-500" />
            3. Partial-Wave Expansions & Interferometry
          </div>
          <p className="text-[#94a3b8]">
            At low impact kinetic energies, de Broglie wavelengths match atomic shell sizes. Quantized electron wavefunctions undergo **constructive and destructive phase interference** inside the atomic potential well.
          </p>
          <p className="text-[#94a3b8]">
            This wave phenomenon creates the characteristic **diffraction peaks and valleys** visible in the plots at larger angles. High-<code className="bg-[#090a0f] px-1 py-0.5 rounded text-[#dfba73] border border-[#222430] font-mono">Z</code> elements produce distinct multiple peaks due to deep potential resonances.
          </p>
        </div>

        {/* Card 4: Sherman Function */}
        <div className="flex flex-col gap-2.5">
          <div className="flex items-center gap-1.5 text-white font-serif font-bold text-[13px] tracking-wide">
            <HelpCircle className="w-4.5 h-4.5 text-emerald-450" />
            4. Sherman Spin S(θ) Polarisation
          </div>
          <p className="text-[#94a3b8]">
            Unpolarized incident beams become **spin-polarized** upon scatter, caused by spin-orbit interactions ( relativistic Mott effects) near the atomic nucleus.
          </p>
          <p className="text-[#94a3b8]">
            The **Sherman function S(θ)** measures this left-right asymmetry, ranging strictly from -1 to +1. High peak Sherman absolute values signify ideal angles for constructing polarized scientific sources!
          </p>
        </div>
      </div>

      <div className="bg-[#090a0f] text-[#94a3b8] p-5 rounded-xl border border-[#222430] mt-3 flex items-start gap-4 select-none leading-normal">
        <div className="font-mono text-center tracking-tighter text-[#c5a059] shrink-0 select-none">
          <div className="text-xl font-bold bg-[#1c1d26] px-3 py-2 rounded-lg border border-[#c5a059]/20 shadow-inner">
            E.S.
          </div>
        </div>
        <div className="text-xs leading-relaxed font-sans">
          <div className="font-bold text-[12.5px] text-white mb-1 font-serif tracking-wide">
            eScatter ELSEPA Computational Lab
          </div>
          Our engine performs a continuous summation of up to 40 Legendre Polynomial orders in real-time, solving the relativistic scattering complex components. This matches the exact analytical and wave-interference boundaries defined inside the official ELSEPA Fortran distribution.
        </div>
      </div>
    </div>
  );
}
