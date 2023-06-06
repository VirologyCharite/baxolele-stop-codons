def adjustInitialReferenceOffsets(offsets):
    result = offsets[:]
    leadingNoneCount = 0
    firstNonNoneOffset = None

    for offset in result:
        if offset is None:
            leadingNoneCount += 1
        else:
            firstNonNoneOffset = offset
            break

    for i in range(leadingNoneCount):
        result[i] = -1 * (firstNonNoneOffset - leadingNoneCount + i)

    return result


def adjustFinalReferenceOffsets(offsets):
    result = offsets[::-1]
    trailingNoneCount = 0
    lastNonNoneOffset = None

    for offset in result:
        if offset is None:
            trailingNoneCount += 1
        else:
            lastNonNoneOffset = offset
            break

    for i in range(trailingNoneCount):
        result[i] = -1 * (lastNonNoneOffset + trailingNoneCount - i)

    return result[::-1]
