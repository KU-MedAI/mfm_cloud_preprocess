"""기본36채널 생성과 기존12채널 검증의 공용 진입점임."""
from pathlib import Path
import argparse
import json
import pysam
from . import MODE
from .schema import CHANNELS, FEATURE_SET, LAYOUT, SCHEMA, LEGACY_SCHEMA
from .complete_export import export, stamp
from .legacy_export import filter_config, reference_contract, hap_reference_contract, _projection
from .variants import VariantLookup


def validate_contig(directory, *, source=False, metrics=None):
    directory = Path(directory)
    metadata = json.loads((directory/'contig_metadata.json').read_text())
    schema = metadata.get('pileup', {}).get('provenance', {}).get('schema')
    if schema == LEGACY_SCHEMA:
        from .legacy_export import validate_contig as validate
    elif schema == SCHEMA:
        from .complete_export import validate_contig as validate
    else:
        raise ValueError('지원하지 않는 pileup schema: '+str(schema))
    return validate(directory, source=source, metrics=metrics)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'alignment', 'reference', 'alignment-reference', 'variants'):
        parser.add_argument('--'+name, required=True, type=Path)
    parser.add_argument('--variant-sample')
    parser.add_argument('--read-type', required=True, choices=['long_read', 'short_read'])
    parser.add_argument('--min-mapq', type=int, default=30)
    parser.add_argument('--contigs', required=True)
    parser.add_argument('--check-inputs', action='store_true')
    args = parser.parse_args()
    contigs = args.contigs.replace(',', ' ').split()
    sample = args.variant_sample or args.root.name
    if args.check_inputs:
        if not any(Path(str(args.variants)+suffix).is_file() for suffix in ('.tbi', '.csi')):
            raise ValueError('BAF sequence VCF index가 필요함')
        lookup = VariantLookup(args.variants, sample)
        lookup.close()
        with pysam.AlignmentFile(str(args.alignment), 'r', reference_filename=str(args.alignment_reference)) as bam:
            if not bam.has_index():
                raise ValueError('reference alignment index가 필요함')
            for contig in contigs:
                contract = reference_contract(bam, args.reference, args.alignment_reference, contig)
                hap_reference_contract(args.root, bam, args.reference, contig)
                for hap in ('HAP1', 'HAP2'):
                    _projection(args.root, contig, hap, contract['length'])
        result = dict(status='pass', scope='reference/chain/VCF identity', contigs=contigs,
                      feature_set=FEATURE_SET, channels=len(CHANNELS))
    else:
        result = export(args.root, args.alignment, args.reference, args.alignment_reference,
                        variants=args.variants, variant_sample=sample,
                        contigs=contigs, read_type=args.read_type, min_mapq=args.min_mapq)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
