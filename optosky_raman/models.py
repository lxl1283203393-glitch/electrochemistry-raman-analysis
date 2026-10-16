"""Immutable provenance with independently retained processing arrays."""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any
import numpy as np


@dataclass
class SpectrumRecord:
    path: Path
    batch: str
    x: np.ndarray
    channels: dict[str, np.ndarray]
    metadata: dict[str, Any] = field(default_factory=dict)
    state: str = 'unknown'
    potential_v: float | None = None
    acquired_at: datetime | None = None
    warnings: list[str] = field(default_factory=list)
    sha256: str = ''
    encoding: str = ''
    stages: dict[str, np.ndarray] = field(default_factory=dict)
    history: list[dict[str, Any]] = field(default_factory=list)

    @property
    def id(self) -> str:
        relative = self.metadata.get('source_relative_path', self.path.name)
        return f'{self.batch}/{relative}'

    def to_ramanspy(self, channel: str = 'dark_subtracted'):
        """Return RamanSPy Spectrum; RamanSPy sorts its own copied axis."""
        import ramanspy as rp
        if channel in self.channels:
            y = self.channels[channel]
        elif channel in self.stages:
            y = self.stages[channel]
        else:
            raise KeyError(f'{self.id}: channel/stage {channel!r} unavailable')
        return rp.Spectrum(np.asarray(y).copy(), np.asarray(self.x).copy())
